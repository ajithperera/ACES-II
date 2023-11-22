#include "flags.h"
      Subroutine Pccd_frmblocks(grad,grad_oo,grad_vv,grad_vo,
     +                          grad_ov,work,Imemleft,Lenoo,
     +                          Lenvv,Lenvo,Nocc,Nvrt,Nbas)

      Implicit Double Precision (A-H, O-Z)

      Dimension Grad(Nbas,Nbas)
      Dimension Grad_oo(Lenoo),Grad_vv(Lenvv)
      Dimension Grad_vo(Lenvo),Grad_ov(Lenvo)
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
      Ovcount = 0
      Length  = Nbas*Nbas

#ifdef _DEBUG_LVL0
      write(6,*)
      write(6,"(a,a)") "@-pccd_frmblocks: The incomming full matrix"
      call output(Grad,1,Nbas,1,Nbas,Nbas,Nbas,1)
#endif

C The occupied/occupied block

      Do Iocc = 1, Nocc
         Do Jocc = 1, Nocc
            Oocount = Oocount + 1
            Grad_oo(Oocount) = Grad(Jocc,Iocc)
         Enddo
      Enddo

C The virtual/virtual block

      Do Ivrt = 1, Nvrt
         Do Jvrt = 1, Nvrt
            Vvcount = Vvcount + 1
            Grad_vv(Vvcount) = Grad(Jvrt+Nocc,Ivrt+Nocc)
         Enddo
      Enddo

C The virtual-occupied (Vo)block

      Do Iocc = 1, Nocc
         Do Jvrt = 1, Nvrt
            Vocount = Vocount + 1
            Grad_vo(Vocount) = Grad(Jvrt+Nocc,Iocc)
         Enddo
      Enddo

C The occupied-virtual (ov) block

      Do Ivrt = 1, Nvrt
         Do Jocc = 1, Nocc
             Ovcount = Ovcount + 1
             Grad_ov(Ovcount) = Grad(Jocc,Ivrt+Nocc)
         Enddo
      Enddo

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
      write(6,*)
      write(6,"(a,a)") "@-pccd_frmblocks: The incomming full matrix"
      call output(Grad,1,Nbfns,1,Nbfns,Nbfns,Nbfns,1)
#endif
C The occupied/occupied block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Iocc = 1, Pop(Irrepr,Ispin)
               Do Jocc = 1, Pop(Irrepl,Ispin)
                  Oocount = Oocount +1
                  Grad_oo(Oocount) = (Grad
     +                                   (Jocc+Ioffo(Irrepl),
     +                                    Iocc+Ioffo(Irrepr)))
               Enddo
            Enddo
         Enddo

C The virtual/virtual block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Ivrt = 1, Vrt(Irrepr,Ispin)
               Do Jvrt = 1,Vrt(Irrepl,Ispin)
                  Vvcount = Vvcount +1
                  Grad_vv(Vvcount) = (Grad
     +                                   (Jvrt+Ioffv(Irrepl),
     +                                    Ivrt+Ioffv(Irrepr)))
               Enddo
            Enddo
         Enddo

C The virtual-occupied (Vo)block

         Do Irrepr = 1, Nirrep
            Irrepl = Dirprd(Irrepr,Irrepx)
            Do Iocc = 1, Pop(Irrepr,Ispin)
               Do Jvrt = 1, Vrt(Irrepl,Ispin)
                  Vocount = Vocount +1
                  Grad_vo(Vocount) = (Grad
     +                                   (Jvrt+IoffV(Irrepl),
     +                                    Iocc+Ioffo(Irrepr)))
               Enddo
            Enddo
         Enddo

C The occupied-virtual (ov) block

         Do Irrepr = 1, Nirrep
           Irrepl = Dirprd(Irrepr,Irrepx)
           Do Jocc = 1, Pop(Irrepl,Ispin)
               Do Ivrt = 1, Vrt(Irrepr,Ispin)
                  Ovcount = Ovcount + 1
                  Grad_ov(Ovcount) = (Grad
     +                                   (Jocc+Ioffo(Irrepl),
     +                                    Ivrt+Ioffv(Irrepr)))
               Enddo
            Enddo
         Enddo

       Enddo 
#endif 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The OO-MO block of the matrix"
      call output(Grad_oo,1,nocc,1,nocc,nocc,nocc,1)
      Write(6,"(a)") "The VV-MO block of the matrix"
      call output(Grad_vv,1,nvrt,1,nvrt,nvrt,nvrt,1)
      Write(6,"(a)") "The VO-MO block of the matrix"
      call output(Grad_vo,1,nvrt,1,nocc,nvrt,nocc,1)
      Write(6,"(a)") "The OV-MO block of the matrix"
      call output(Grad_ov,1,nocc,1,nvrt,nocc,nvrt,1)
#endif
      Return
      End
