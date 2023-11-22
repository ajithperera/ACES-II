#include "flags.h"
      Subroutine Tdcc_form_mu_dot_0(Work,Memleft,Irrepx,Iuhf,Mu0_d)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Double Precision Mu0_d
      Integer Hbr_hp_list,Mu0_hp_list

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
       
      Hbr_hp_list = 93
      Mu0_hp_list = 390
      
      Max_vo = Max(Irpdpd(Irrepx,9),Irpdpd(Irrepx,10))

      I000 = 1
      I010 = I000 + Max_vo
      Iend = I010 + Max_vo 
      If (Iend .Gt. Memleft) Call Insmem("@-Tdvee_form_mu_0_0",
     +                                    Memleft,Iend)

C NDeP Eqn. 24, first-term, m(a,i)*Hbar(a,i)

      Do Ispin = 1, 1+Iuhf 
       
         Mu0_d = 0.0D0
         Lenvo = Irpdpd(Irrepx,8+Ispin)

         Call Getlst(Work(I000),1,1,1,Ispin,Mu0_hp_list)
         Call Getlst(Work(I010),1,1,1,Ispin,Hbr_hp_list)

         Mu0_d =  Mu0_d + Ddot(Lenvo,Work(I000),1,Work(1010),1)

#ifdef _DEBUG_LVL0
        Write(6,*) 
CSSS        call checksum("Tdcc_form_mu_dot_0:",Work(I000),Lenvo,s)
CSSS        call checksum("Tdcc_form_mu_dot_0:",Work(I010),Lenvo,s)
        if (Ispin .Eq.1) Write(6,"(a,F12.6)") "The mu_0_dot alpha:",
     +                   Mu_0_d
        if (Ispin .Eq.2) Write(6,"(a,F12.6)") "The mu_0_dot  beta:",
     +                   Mu_0_d
#endif

      Enddo 

      Return
      End
