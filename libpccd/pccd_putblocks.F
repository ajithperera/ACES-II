#include "flags.h"
      Subroutine Pccd_Putblocks(Grad,Grad_oo,Grad_vv,Grad_vo,
     +                          Grad_ov,work,Imemleft,Lenoo,
     +                          Lenvv,Lenvo,Nocc,Nvrt,Nbas,
     +                          String,Sym_packed)

      Implicit Double Precision (A-H, O-Z)
      Character*7 String 
      Logical Symmetry,Sym_packed 

      Dimension Grad(Nbas,Nbas)
      Dimension Grad_oo(Lenoo),Grad_vv(Lenvv)
      Dimension Grad_vo(Lenvo),Grad_ov(Lenvo)
      Dimension Work(Imemleft)
      Dimension Ioffo(8),Ioffv(8)
      Integer Oocount,Vvcount,Vocount,Ovcount

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "info.com"
#include "sym.com"

      Common/Symm/Symmetry

      Data Ione /1/

      Irrepx  = 1
      Oocount = 0
      Vvcount = 0
      Vocount = 0
      Ovcount = 0
      
      Call Getrec(0,"JOBARC","SPNINDEX",Len,Ispin)
      If (Len .Gt .0) Then 
         Call Getrec(20,"JOBARC","SPNINDEX",Ione,Ispin)
      Else
         Ispin = Ione
      Endif 

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
      write(6,*)
      write(6,"(a,i2)") "@-pccd_putblocks, The Ispin:", Ispin
      write(6,"(a,1x,8i2)") "ioffo:",(Ioffo(i),i=1,Nirrep)
      write(6,"(a,1x,8i2)") "ioffv:",(Ioffv(i),i=1,Nirrep)

      Write(6,*)
      Write(6,"(a)") "Printing from Pccd_putblocks"
      Write(6,"(a)") "The OO-MO gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_oo(i),i=1,Lenoo)
      Write(6,"(a)") "The VV-MO  gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_vv(i),i=1,Lenvv)
      Write(6,"(a)") "The VO-MO  gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_vo(i),i=1,Lenvo)
      Write(6,"(a)") "The OV-MO  gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_ov(i),i=1,Lenvo)
#endif

      If (Sym_packed) Then
          Ioffoo = Ione
          Ioffvv = Ione
          Ioffvo = Ione

          Nocc2 = Nocc*Nocc
          Nvrt2 = Nvrt*Nvrt
          Nvovo = Nvrt*Nocc
          Call Dzero(Grad,Nbas*Nbas)

          I000 = Ione
          I010 = I000 + Nocc2 
          Iend = I010 + Nocc2 
          Call Putblk(Grad,Grad_oo(Ioffoo),"OO",Work(I000),
     +                Work(I010),Nbas,Irpdpd(Irrepx,20+Ispin),
     +                Nocc2,Ispin)
          I000 = Ione
          I010 = I000 + Nvrt2 
          Iend = I010 + Nvrt2
          Call Putblk(Grad,Grad_vv(Ioffvv),"VV",Work(I000),
     +                Work(I010),Nbas,Irpdpd(Irrepx,18+Ispin),
     +                Nvrt2,Ispin)
          I000 = Ione
          I010 = I000 + Nvovo
          Iend = I010 + Nvovo
          Call Putblk(Grad,Grad_vo(Ioffvo),"VO",Work(I000),
     +                Work(I010),Nbas,Irpdpd(Irrepx,8+Ispin),
     +                Nvovo,Ispin)

          If (String .Eq. "Special") Then

          Do Irrepr = 1, Nirrep
             Irrepl = Dirprd(Irrepr,Irrepx)
             Do Iocc = 1, Pop(Irrepr,Ispin)
                Do Ivrt = 1, Vrt(Irrepl,Ispin)
                  Ovcount =  Ovcount + 1 
                  Grad(iocc+ioffo(irrepr),ivrt+ioffv(irrepl)) = 
     +                                    Grad_ov(Ovcount)
                Enddo
             Enddo
          Enddo

          Else

          I000 = Ione
          I010 = I000 + Nvovo
          Iend = I010 + Nvovo
          Call Putblk(Grad,Grad_ov(Ioffvo),"OV",Work(I000),
     +                Work(I010),Nbas,Irpdpd(Irrepx,8+Ispin),
     +                Nvovo,Ispin)
          Endif 

          Ioffo  = Ioffo  + Irpdpd(Irrepx,20+Ispin)
          Ioffv  = Ioffv  + Irpdpd(Irrepx,18+Ispin)
          Ioffvo = Ioffvo + Irpdpd(Irrepx,8+Ispin)
             
      Else

C The occupied/occupied block

          Call Dzero(Grad,Nbas*Nbas)
          Do Iocc = 1, Nocc
             Do Jocc = 1, Nocc
                 Oocount = Oocount +1
                 Grad(Jocc,Iocc) = Grad_oo(Oocount)
             Enddo
          Enddo 

C The virtual/virtual block

          Do Ivrt = 1, Nvrt
             Do Jvrt = 1,Nvrt
                Vvcount = Vvcount +1
                Grad(Jvrt+Nocc,Ivrt+Nocc) = Grad_vv(Vvcount)
            Enddo
          Enddo 

C The virtual-occupied (Vo)block

          Do Iocc = 1, Nocc
             Do Jvrt = 1, Nvrt
                Vocount = Vocount + 1
                Grad(Jvrt+Nocc,Iocc) = Grad_vo(Vocount)
             Enddo
          Enddo

C The occupied-virtual (ov) block

          If (String .Eq. "Special") Then

             Do Jocc = 1, Nocc
                Do Ivrt = 1, Nvrt
                   Ovcount = Ovcount + 1
                   Grad(Jocc,Ivrt+Nocc) = Grad_ov(Ovcount)
                Enddo
             Enddo

          Else

             Do Ivrt = 1, Nvrt
                Do Jocc = 1, Nocc
                   Ovcount = Ovcount + 1
                   Grad(Jocc,Ivrt+Nocc) = Grad_ov(Ovcount)
                Enddo
             Enddo 

          Endif 
      ENdif 

#ifdef _NOSKIP

          Do Ispin = Ione, Ione 

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
      Write(6,*)
      Write(6,"(a)") "Printing from Pccd_putblocks"
      Write(6,"(a)") "The OO-MO antisym. gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_oo(Ioff_oo-1+i),
     +                          i=1,Ndim_oo)
      Write(6,"(a)") "The VV-MO antisym. gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_vv(Ioff_vv-1+i),
     +                          i=1,Ndim_vv)
      Write(6,"(a)") "The VO-MO antisym. gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_vo(Ioff_vo-1+i),
     +                          i=1,Ndim_vo)
      Write(6,"(a)") "The OV-MO antisym. gradient matrices"
      write(6,"(6(1x,F12.6))") (Grad_ov(Ioff_vo-1+i),
     +                          i=1,Ndim_vo)
#endif
             Call Dzero(Grad,Nbas*Nbas)

C The occupied/occupied block

             Do Irrepr = 1, Nirrep
                Irrepl = Dirprd(Irrepr,Irrepx)
                Do Iocc = 1, Pop(Irrepr,Ispin)
                   Do Jocc = 1, Pop(Irrepl,Ispin)
                      Oocount = Oocount +1
                      Grad(Jocc+Ioffo(Irrepl),Iocc+Ioffo(Irrepr)) 
     +                                        = Grad_oo(Oocount)
                   Enddo
                Enddo
             Enddo

C The virtual/virtual block

            Do Irrepr = 1, Nirrep
               Irrepl = Dirprd(Irrepr,Irrepx)
               Do Ivrt = 1, Vrt(Irrepr,Ispin)
                  Do Jvrt = 1,Vrt(Irrepl,Ispin)
                     Vvcount = Vvcount +1
                     Grad(Jvrt+Ioffv(Irrepl),Ivrt+Ioffv(Irrepr))
     +                                      = Grad_vv(Vvcount)
                  Enddo
               Enddo
            Enddo

C The virtual-occupied (Vo)block

            Do Irrepr = 1, Nirrep
               Irrepl = Dirprd(Irrepr,Irrepx)
               Do Iocc = 1, Pop(Irrepr,Ispin)
                  Do Jvrt = 1, Vrt(Irrepl,Ispin)
                     Vocount = Vocount +1
                    Grad(Jvrt+IoffV(Irrepl),Iocc+Ioffo(Irrepr))
     +                                      = Grad_vo(Vocount)
                  Enddo
               Enddo
            Enddo

C The occupied-virtual (ov) block

            If (String .Eq. "Ov_like") Then 
               Do Irrepr = 1, Nirrep
                  Irrepl = Dirprd(Irrepr,Irrepx)
                  Do Ivrt = 1, Vrt(Irrepr,Ispin)
                     Do Jocc = 1, Pop(Irrepl,Ispin)
                        Ovcount = Ovcount + 1
                        Grad(Jocc+Ioffo(Irrepl),Ivrt+Ioffv(Irrepr))
     +                                          = Grad_ov(Ovcount)
                     Enddo
                  Enddo
               Enddo
             Endif 

            If (String .Eq. "Vo_like") Then
                Do Irrepr = 1, Nirrep
                   Irrepl = Dirprd(Irrepr,Irrepx)
                   Do Jocc = 1, Pop(Irrepl,Ispin)
                      Do Ivrt = 1, Vrt(Irrepr,Ispin)
                         Ovcount = Ovcount + 1
                        Grad(Jocc+Ioffo(Irrepl),Ivrt+Ioffv(Irrepr))
     +                                          = Grad_ov(Ovcount)
                      Enddo
                   Enddo
                Enddo
            Endif

          Enddo 
      Endif 
#endif 

      Return
      End
