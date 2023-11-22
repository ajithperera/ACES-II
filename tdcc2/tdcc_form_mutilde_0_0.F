#include "flags.h"
      Subroutine Tdcc_form_mutilde_0_0(Work,Memleft,Doo,Dvv,Dvo,
     +                                 Irrepx,Lenoo,Lenvv,Lenvo,
     +                                 Iuhf,Mu_t_0)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Dimension Doo(Lenoo),Dvv(Lenvv),Dvo(Lenvo)
      Double Precision Mu_ii_0,Mu_t_l,Mu_t_0

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "sym.com"

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a)") "---Entered tdcc_form_mutilde_0_0---"
        Write(6,*)
#endif
      I000 = 1
      Iend = I000 + Lenvo 
      If (Iend .Gt. Memleft) Call Insmem("@-Tdvee_form_mu_0_0",
     +                                    Memleft,Iend)

C Both contributions contribute only when irrepx=1.
C NDeP Eqn. 19, first-term, M(I,J)*delta(I,J) + M(i,j)*delta(i,j)

      L1_list = 190 
      Ioff_vo = 1

      Do Ispin = 1, 1+Iuhf 

         Mu_ii_0 = 0.0D0
         Mu_t_l  = 0.0D0
         Mu_t_0  = 0.0D0

         Istart = (Ispin-1)*Irpdpd(Irrepx,21)

         Do Irrep_j = 1, Nirrep
            Irrep_i = Dirprd(Irrep_j,Irrepx)
            If (Irrep_i .EQ. Irrep_j) Then
               Do Jocc = 1, Pop(Irrep_j,Ispin)
                  Ioff = Istart + Jocc +  
     +                   (Jocc-1)* Pop(Irrep_i,Ispin)
                   Mu_ii_0 = Mu_ii_0 + Doo(Ioff) 
               Enddo
            Endif
            Istart = Ioff
         Enddo

C NDeP Eqn. 19, second-term, Mu(A,I)*L(A,I) + Mu(a,i)*L(a,i)

         If (Irrepx .EQ. 1) Then
            Ioff_vo = Ioff_vo + (Ispin-1) * Irpdpd(Irrepx,9)
            Ndim_vo =  Irpdpd(Irrepx,8+ispin)

            Call Getlst(Work(I000),1,1,1,Ispin,L1_list)
            Mu_t_l =  Ddot(Ndim_vo,Work(I000),1,Dvo(ioff_vo),1)
            Mu_t_0 = Mu_ii_0 + Mu_t_l

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") "L1(i,a) at list 190"
      Write(6,"(6(1x,F12.6))") (Work(I000+i-1),i=1,Ndim_vo)
      Write(6,"(a)") "Mu(a,i)"
      Write(6,"(6(1x,F12.6))") (Dvo(ioff_vo+i-1),i=1,Ndim_vo)
#endif 
         Endif 

#ifdef _DEBUG_LVL0
      If (Ispin .Eq.1) Write(6,"(a,a,3(1x,F12.8))") "The mu_ii_0,",
     +                 "mu_t_l,mu_t_0, alpha", Mu_ii_0,Mu_t_l,
     +                  Mu_t_0
      If (Ispin .Eq.2) Write(6,"(a,a,3(1x,F12.8))") "The mu_ii_0,",
     +                 "mu_t_l,mu_t_0, beta ", Mu_ii_0,Mu_t_l,
     +                  Mu_t_0
#endif
      Enddo 

      Return
      End
