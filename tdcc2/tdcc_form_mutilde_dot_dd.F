      Subroutine Tdcc_form_mutilde_dot_dd(Work,Memleft,Irrepx,Iuhf,
     +                               Mut_dot)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)
      Double Precision Mut_dot_ab, Mut_dot_pp
      Integer D2_list
      Logical Rhf

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"

      Rhf = (Iuhf.Eq.0)
      Mut_dot_ab = 0.0D0

C Mut_dot = Mut(Ab,Ij)*(e_I+e_j-e_A-e_b)

      Mut_list = 326
      D2_list  = 66

      Do Irrep_ij=1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)

c Since Irrepx can only be one; irrep_ab=Irrep_ij

            Ndim_ij_z = Irpdpd(Irrep_ij,14)
            Ndim_ab_z = Irpdpd(Irrep_ab,13)
            Ndim_ij_d = Ndim_ij_z
           
            Call Getlst(Work(I000),1,Ndim_ij_z,1,Irrep_ij,Mut_list)
            Call Getlst(Work(I010),1,Ndim_ij_d,1,Irrep_ij,D2_list)
        
            Mut_dot_ab = Mut_dot_ab + Ddot(Ndim_ij_z*Ndim_ab_z,
     +                                Work(I000),1,Work(I010),1)
  
      Enddo

C Mut_dot = Mut(AB,IJ)*(e_I+e_J-e_A-e_B) ISPIN=1
C Mut_dot = Mut(ab,ij)*(e_i+e_j-e_a-e_b) ISPIN=2
C Read Mut as Mut(A<B,I<J) and D(A<B,I<J)

#ifdef _DEBUG_LVL0
      If (Ispin .Eq. 1) Write(6,"(a,1x,F12.6)")
     +                   "@-Tdcc_form_mutilde_dot_dd, AbAb:",
     +                       Mut_dot_ab
#endif 
      If (Rhf) Return 
   
      Mut_dot_aa = 0.0D0
      Mut_dot_bb = 0.0D0

      Do Ispin = 1, 1+Iuhf

         Mut_list = 323 + Ispin 
         D2_list  = 63  + Ispin
 
         Do Irrep_ij=1, Nirrep
            Irrep_ab = Dirprd(Irrep_ij,Irrepx)

c Since Irrepx can only be one; irrep_ab=Irrep_ij

            Ndim_ij_z = Irpdpd(Irrep_ij,2+Ispin)
            Ndim_ab_z = Irpdpd(Irrep_ab,Ispin)
            Ndim_ij_d = Ndim_ij_z
         
            Call Getlst(Work(I000),1,Ndim_ij_z,1,Irrep_ij,Mut_list)
            Call Getlst(Work(I010),1,Ndim_ij_d,1,Irrep_ij,D2_list)
#ifdef _DEBUG_LVLM
           call output(Work(I010),1,Ndim_ab_z,1,Ndim_ij_z,Ndim_ab_z,
     +                  Ndim_ij_z,1)
#endif 

            If (Ispin .EQ. 1) Then
               Mut_dot_aa = Mut_dot_aa + Ddot(Ndim_ij_z*Ndim_ab_z,
     +                                   Work(I000),1,Work(I010),1)
            Else
               Mut_dot_bb = Mut_dot_bb + Ddot(Ndim_ij_z*Ndim_ab_z,
     +                                   Work(I000),1,Work(I010),1)
            Endif 
         Enddo 
         
#ifdef _DEBUG_LVL0
      Write(6,*)
      If (Ispin .Eq. 1) Write(6,"(a,1x,F12.6)")
     +                   "@-Tdcc_form_mutilde_dot_dd, AAAA:",
     +                    Mut_dot_aa
      If (Ispin .Eq. 1) Write(6,"(a,1x,F12.6)")
     +                   "@-Tdcc_form_mutilde_dot_dd, BBBB:",
     +                    Mut_dot_bb
#endif 

      Enddo 
      Return
      End
    
