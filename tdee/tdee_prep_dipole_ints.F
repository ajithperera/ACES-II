#include "flags.h"
      Subroutine Tdee_prep_dipole_ints(Dip_ints,Dip_ints_full,Scr1,
     +                              Scr2,Dip_ints_oo,Dip_ints_vv,
     +                              Dip_ints_vo,Work,Imemleft,Lenoo,
     +                              Lenvv,Lenvo,Nbfns,Naobfns,Iuhf,
     +                              Ipert,Label,Irrepx)

      Implicit Double Precision (A-H, O-Z)

      Dimension Scr1(Naobfns*Nbfns),Scr2(Naobfns*Nbfns)
      Dimension Dip_ints(Naobfns*(Naobfns+1)/2)
      Dimension Dip_ints_full(Nbfns,Nbfns)
      Dimension Dip_ints_oo(Lenoo),Dip_ints_vv(Lenvv)
      Dimension Dip_ints_vo(Lenvo)
      Dimension Work(Imemleft)
      Dimension Ioffo(8),Ioffv(8)
      Character*8 Label 
      Integer Oocount,Vvcount,Vocount

      Logical Vprop, Abacus 

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"
#include "info.com"
#include "sym.com"

      Vprop   = .True.
      Abacus  = .False.
      Length  = (Naobfns+1)*Naobfns/2
      Oocount = 0
      Vvcount = 0
      Vocount = 0

      Do Ispin = 1, Iuhf+1

         Ioff_oo = 1 + (Ispin-1) * Irpdpd(Irrepx,21)
         Ioff_vv = 1 + (Ispin-1) * Irpdpd(Irrepx,19)
         Ioff_vo = 1 + (Ispin-1) * Irpdpd(Irrepx,9)

         Ndim_oo = Irpdpd(Irrepx,20+Ispin)
         Ndim_vv = Irpdpd(Irrepx,18+Ispin)
         Ndim_vo = Irpdpd(Irrepx,8+Ispin)

C Calculate offsets occ and vrts.

         Ioffo(1) = 0
         Ioffv(1) = Nocco(Ispin)
      
         Do Irrep = 1, Nirrep-1
            Ioffo(Irrep+1) = Ioffo(Irrep) + Pop(Irrep,Ispin)
            Ioffv(Irrep+1) = Ioffv(Irrep) + Vrt(Irrep,Ispin)
         Enddo

         If (Vprop) Then 
            Call Getrec(20,"JOBARC",Label,Length*Iintfp,
     +                 Dip_ints)
            Call Expnd2(Dip_ints,Dip_ints_full,Naobfns)
            Call Ao2mo2(Dip_ints_full,Dip_ints_full,Scr1,Scr2,
     +                    Nbfns,Naobfns,Ispin)

         Else if (Abacus) Then

            Call Getlst(Dip_ints,Ipert,1,1,Irrepx, 198)
            Call Expnd2(Dip_ints,Dip_ints_full,Naobfns)
            Call Ao2mo2(Dip_ints_full,Dip_ints_full,Scr1,Scr2,
     +                    Nbfns,Naobfns,Ispin)
         Endif
#ifdef _DEBUG_LVL0
      write(6,"(a)") "The MO basis dipole integrals"
      call output(Dip_ints_full,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
#endif

C The occupied/occupied block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Iocc = 1, Pop(Irrepr,Ispin)
               Do Jocc = 1, Pop(Irrepl,Ispin)
                  Oocount = Oocount +1
                  Dip_ints_oo(Oocount) = Dip_ints_full
     +                                           (Jocc+Ioffo(Irrepl),
     +                                            Iocc+Ioffo(Irrepr))
               Enddo
            Enddo
         Enddo 

C The virtual/virtual block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Ivrt = 1, Vrt(Irrepr,Ispin)
               Do Jvrt = 1,Vrt(Irrepl,Ispin)
                  Vvcount = Vvcount +1
                  Dip_ints_vv(Vvcount) = Dip_ints_full
     +                                           (Jvrt+Ioffv(Irrepl),
     +                                            Ivrt+Ioffv(Irrepr))
               Enddo
            Enddo
         Enddo 

C The virtual/occupied block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Iocc = 1, Pop(Irrepr,Ispin)
               Do Jvrt = 1, Vrt(Irrepl,Ispin)
                  Vocount = Vocount +1
                  Dip_ints_vo(Vocount) = Dip_ints_full
     +                                           (Jvrt+IoffV(Irrepl),
     +                                            Iocc+Ioffo(Irrepr))
               Enddo
            Enddo
         Enddo 

C Save the MO basis,sym-packed dipole integral blocks 

         Call Putlst(Dip_ints_oo(Ioff_oo),1,1,1,Ispin,370)
         Call Putlst(Dip_ints_vv(Ioff_vv),1,1,1,Ispin,371)
         Call Putlst(Dip_ints_vo(Ioff_vo),1,1,1,Ispin,372)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The OO-MO dipole matrices"
      write(6,"(6(1x,F12.6))") (Dip_ints_oo(Ioff_oo-1+i),
     +                          i=1,Ndim_oo)
      Write(6,"(a)") "The VV-MO dipole matrices"
      write(6,"(6(1x,F12.6))") (Dip_ints_vv(Ioff_vv-1+i),
     +                          i=1,Ndim_vv)
      Write(6,"(a)") "The VO-MO dipole matrices"
      write(6,"(6(1x,F12.6))") (Dip_ints_vo(Ioff_vo-1+i),
     +                          i=1,Ndim_vo)
#endif 

      Enddo
         
      Return
      End
