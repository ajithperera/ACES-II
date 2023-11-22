#include "flags.h"
      Subroutine Tdcc_form_mu_0_0(Work,Memleft,Doo,Dvv,Dvo,
     +                            Irrepx,Lenoo,Lenvv,Lenvo,
     +                            Iuhf,Mu_0)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)
      Double Precision Mu_0

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "sym.com"

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered tdcc_form_mu_0_0---"
        Write(6,*)
#endif
      I000 = 1
      Iend = I000 + Lenvo 
      If (Iend .Gt. Memleft) Call Insmem("@-Tdvee_form_mt_0_0",
     +                                    Memleft,Iend)
#ifdef _DEBUG_LVLM
      do ispin=1, 1+iuhf
      ioff = 1 + (Ispin-1)*Irpdpd(Irrepx,21)
      do irrep_i=1,Nirrep
         irrep_j=dirprd(irrep_i,Irrepx)
         ndim_i = Pop(irrep_i,ispin)
         Ndim_j = Pop(irrep_j,ispin)
         write(6,*) irrep_i, Irrep_j,ioff,Ndim_i,Ndim_j
         call output(doo(ioff),1,ndim_j,1,Ndim_i,Ndim_j,Ndim_i,1)
         Ioff = Ioff + ndim_i*Ndim_j
      enddo
      enddo
#endif 
C This can contribute only when irrepx=1
C Mu(i,j)Delta(i,j) + Mu(I,J)Delta(I,J)

      Do Ispin = 1, 1+Iuhf 
         Mu_0 = 0.0D0
      
         Istart = (Ispin-1)*Irpdpd(Irrepx,21)

         Do Irrep_j = 1, Nirrep
            Irrep_i = Dirprd(Irrep_j,Irrepx)
            If (Irrep_i .EQ. Irrep_j) Then
                Do Jocc = 1, Pop(Irrep_j,Ispin)
                   Ioff = Istart + Jocc +  
     +                    (Jocc-1)*Pop(Irrep_i,Ispin)
#ifdef _DEBUG_LVLM
               Write(6,"(3a,6(1x,I4))") "irrep_i,Irrep_j,",
     +                   "Pop(Irrep_j,Ispin),",
     +                   "Pop(Irrep_i,Ispin),istart,Ioff",
     +                   irrep_i,Irrep_j,
     +                   Pop(Irrep_j,Ispin),
     +                   Pop(Irrep_i,Ispin),istart,Ioff
#endif 
                   Mu_0 = Mu_0 + Doo(Ioff) 
                Enddo
            Endif
            Istart = Ioff 
         Enddo

#ifdef _DEBUG_LVL0
        Write(6,*)
        if (Ispin .Eq.1) Write(6,"(a,F12.6)") "The mu_0 alpha", 
     +                   Mu_0 
        if (Ispin .Eq.2) Write(6,"(a,F12.6)") "The mu_0  beta", 
     +                   Mu_0 
#endif 
      Enddo 

      Return
      End
     
       
