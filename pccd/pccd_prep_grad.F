#include "flags.h"
      Subroutine Pccd_prep_grad(grad,grad_oo,grad_vv,grad_vo,
     +                          grad_ov,work,Imemleft,Lenoo,
     +                          Lenvv,Lenvo,Nocc,Nvrt,Nbfns)

      Implicit Double Precision (A-H, O-Z)

      Dimension Grad(Nbfns,Nbfns)
      Dimension Grad_oo(Nocc,Nocc),Grad_vv(Nvrt,Nvrt)
      Dimension Grad_vo(Nvrt,Nocc),Grad_ov(Nocc,Nvrt)
      Dimension Work(Imemleft)
      Dimension Ioffo(8),Ioffv(8)
      Integer Oocount,Vvcount,Vocount

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "sym.com"

      Irrepx  = 1
      Oocount = 0
      Vvcount = 0
      Vocount = 0
      oVcount = 0
      Length  = Nbfns*Nbfns

#ifdef _NOSKIP
      Do Ispin = 1, 1

         Ioff_oo = 1 + (Ispin-1) * Irpdpd(Irrepx,21)
         Ioff_vv = 1 + (Ispin-1) * Irpdpd(Irrepx,19)
         Ioff_vo = 1 + (Ispin-1) * Irpdpd(Irrepx,9)

         Ndim_oo = Irpdpd(Irrepx,20+Ispin)
         Ndim_vv = Irpdpd(Irrepx,18+Ispin)
         Ndim_vo = Irpdpd(Irrepx,8+Ispin)
         Ndim_ov = Ndim_vo 

C Calculate offsets occ and vrts.

         Ioffo(1) = 0
         Ioffv(1) = Nocco(Ispin)
      
         Do Irrep = 1, Nirrep-1
            Ioffo(Irrep+1) = Ioffo(Irrep) + Pop(Irrep,Ispin)
            Ioffv(Irrep+1) = Ioffv(Irrep) + Vrt(Irrep,Ispin)
         Enddo

#ifdef _DEBUG_LVL0
      write(6,"(a)") "The orbital rotation gradients"
      call output(Grad,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
#endif
C The occupied/occupied block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Iocc = 1, Pop(Irrepr,Ispin)
               Do Jocc = 1, Pop(Irrepl,Ispin)
                  Oocount = Oocount +1
                  Grad_oo(Oocount) = Grad
     +                                   (Jocc+Ioffo(Irrepl),
     +                                    Iocc+Ioffo(Irrepr))
               Enddo
            Enddo
         Enddo

C The virtual/virtual block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Ivrt = 1, Vrt(Irrepr,Ispin)
               Do Jvrt = 1,Vrt(Irrepl,Ispin)
                  Vvcount = Vvcount +1
                  Grad_vv(Vvcount) = Grad
     +                                   (Jvrt+Ioffv(Irrepl),
     +                                    Ivrt+Ioffv(Irrepr))
               Enddo
            Enddo
         Enddo

C The virtual/occupied block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Iocc = 1, Pop(Irrepr,Ispin)
               Do Jvrt = 1, Vrt(Irrepl,Ispin)
                  Vocount = Vocount +1
                  Grad_vo(Vocount) = Grad
     +                                   (Jvrt+IoffV(Irrepl),
     +                                    Iocc+Ioffo(Irrepr))
               Enddo
            Enddo
         Enddo

C The occupied/virtual block

         Do Irrepr = 1, Nirrep
           Irrepl = Dirprd(Irrepr,Irrepx)
            Do Ivrt = 1, Vrt(Irrepr,Ispin)
               Do Jocc = 1, Pop(Irrepl,Ispin)
                  Ovcount = Ovcount + 1
                  Grad_ov(Ovcount) = Grad
     +                                   (Jocc+Ioffo(Irrepr),
     +                                    Ivrt+Ioffv(Irrepl))
               Enddo
            Enddo
         Enddo

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The OO-MO gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_oo(Ioff_oo-1+i),
     +                          i=1,Ndim_oo)
      Write(6,"(a)") "The VV-MO gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_vv(Ioff_vv-1+i),
     +                          i=1,Ndim_vv)
      Write(6,"(a)") "The VO-MO gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_vo(Ioff_vo-1+i),
     +                          i=1,Ndim_vo)
      Write(6,"(a)") "The OV-MO gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_ov(Ioff_vo-1+i),
     +                          i=1,Ndim_vo)
#endif
      Enddo
#endif 

      Call Zero(Grad_oo,Nocc*Nocc)
      Call Zero(Grad_vv,Nvrt*Nvrt)
      Call Zero(Grad_vo,Nvrt*Nocc)
      Call Zero(Grad_ov,Nocc*Nvrt)

      Do I = 1, Nocc
         Do J = 1, Nocc
            Grad_oo(J,I) = Grad(J,I)
         Enddo
      Enddo 

      Ioff = 0
      Do I = Nocc+1, Nbfns 
         Ioff = Ioff + 1
         Joff = 0
         Do J = Nocc+1, Nbfns 
            Joff = Joff + 1
            Grad_vv(Joff,Ioff) = Grad(J,I)
         Enddo
      Enddo 

      Ioff = 0
      Do I = 1, Nocc
         Ioff = Ioff + 1
         Joff = 0
         Do J = Nocc+1, Nbfns
            Joff = Joff + 1
            Grad_vo(Joff,Ioff) = Grad(J,I)
         Enddo
      Enddo 

      Ioff = 0
      Do I = 1, Nocc 
         Ioff = Ioff + 1
         Joff = 0
         Do J = Nocc+1, Nbfns 
            Joff = Joff + 1
            Grad_ov(Ioff,Joff) = Grad(I,J)
         Enddo
      Enddo 
#
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The OO-MO gradient matrices"
      call output(Grad_oo,1,Nocc,1,Nocc,Nocc,Nocc,1)
      Write(6,"(a)") "The VV-MO gradient matrices"
      call output(Grad_vv,1,Nvrt,1,Nvrt,Nvrt,Nvrt,1)
      Write(6,"(a)") "The VO-MO gradient matrices"
      call output(Grad_vo,1,Nvrt,1,Nocc,Nvrt,Nocc,1)
      Write(6,"(a)") "The OV-MO gradient matrices"
      call output(Grad_ov,1,Nocc,1,Nvrt,Nocc,Nvrt,1)
#endif 
      Return
      End
