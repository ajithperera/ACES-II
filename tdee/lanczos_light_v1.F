      Subroutine Lanczos_light_v1(Nsize,Vecs,Info)

      Implicit Double Precision (A-H,O-Z)
      Integer Retlabel, Revcom,Info(4)
      Double Precision Ps_mod 

      Parameter (Threshold = 1.0D-10)

      Dimension Vecs(Nsize,6)

#include "lanczos_vars.com"

      Common /Beta/ Beta_jrl,Gama_jrl,J
      Common /Init_norms/Sqrt_Norm,Rl_norm
      Data Retlabel /1/
      Save /lanczos_vars/

      Revcom  = Info(2)
      Info(2) = 0
      If (Revcom .EQ. 0) J = 1

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "--- Entered Lanczos_light ---"
      Write(6,"(a,i4,1x,i2)") "The counter,retlabel :", j,retlabel
#endif
c Save the tridiagoal matrix (Trd defined in lanczos_vars.com) while 
C lanczos iterations.

      Goto (10,30,40) Retlabel 
      
 10   Continue 

      RR_Norm = Ddot(Nsize,Vecs(1,1),1,Vecs(1,1),1)
#ifdef _DEBUG_LVL0
      Write(6,"(a,F15.8)")  "@-Entry; R_Norm : ", Dsqrt(RR_Norm)
#endif
      If (RR_Norm .NE. 0.0D0) Then
        Sqrt_Norm = 1.0D0/Dsqrt(RR_norm)
        Call Dscal(Nsize,Sqrt_Norm,Vecs(1,1),1)
      Else
         Write(6,"(a)") " Division by zero @-lanczos_light_v1"
         call Errex
      Endif

      RL_Norm = Ddot(Nsize,Vecs(1,1),1,Vecs(1,2),1) 
#ifdef _DEBUG_LVL0
      Write(6,"(a,F15.8)")  "@-Entry; RL_Norm: ", RL_Norm
#endif
      If (RL_Norm .NE. 0.0D0) Then
         RL_NORM = 1.0D0/RL_NORM
         Call Dscal(Nsize,RL_NORM,Vecs(1,2),1)
      Else 
         Write(6,"(a)") " Division by zero @-lanczos_light_v1"
         call Errex
      Endif 

C Check whether <R|L> = 1.0 for the starting vectors
     
      RL_Norm = Ddot(Nsize,Vecs(1,1),1,Vecs(1,2),1) 
      If (Dabs(RL_Norm-1.0D0) .Gt. Threshold) Then
         Write(6,"(a)") "The starting vectors are not biorthogonal"
         Call Errex
      Endif 

      Call Dcopy(Nsize,Vecs(1,1),1,Vecs(1,3),1)
      Call Dcopy(Nsize,Vecs(1,2),1,Vecs(1,4),1)

  20  Continue 

C Do |v_j+1> = A |v_j>

      Info(2)  = 1
      Info(3)  = 3
      Info(4)  = 5
      Retlabel = 2 
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "Returning 1"
#endif 
      Return

 30   Continue 

C Do <w_j+1| = <w_j| A

      Info(2)  = 2
      Info(3)  = 4
      Info(4)  = 6
      Retlabel = 3
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "Returning 2"
#endif 
      Return

C Vec(1,1) and Vecs(1,2) keeps v_(j_1) and w_(j_1)
C Vec(1,3) and Vecs(1,4) keeps v_j and w_j
C Vec(1,5) and Vecs(1,6) keeps v_j+1 and w_j+1
      
 40   Continue 

      If (J .eq. 1) Then
      
          Alpa_jrl = -Ddot(Nsize,Vecs(1,5),1,Vecs(1,4),1)

#ifdef _DEBUG_LVL0
          Write(6,"(a,F15.8)") "<L_old|R_new> = ", -Alpa_jrl
#endif

C p_(1) = V_(2) - alpa_(j)V_(1)
      
          Call Daxpy(Nsize,Alpa_jrl,Vecs(1,3),1,Vecs(1,5),1)

C s_(1) = w_(2) - alpa_(j)w_(j)

          Call Daxpy(Nsize,Alpa_jrl,Vecs(1,4),1,Vecs(1,6),1)

          Ps_mod   = Ddot(Nsize,Vecs(1,5),1,Vecs(1,6),1)

          Beta_jp1 = Dsqrt(Dabs(Ps_mod))
          Gama_jp1 = (Ps_mod/Dabs(Ps_mod)) * Beta_jp1

          Beta_jp1_inv = 1.0D0/Beta_jp1
          Gama_jp1_inv = 1.0D0/Gama_jp1
      
          Call Dscal(Nsize,Beta_jp1_inv,Vecs(1,5),1)
          Call Dscal(Nsize,Gama_jp1_inv,Vecs(1,6),1)

	  Beta_jrl = -Beta_jp1
          Gama_jrl = -Gama_jp1

#ifdef _DEBUG_LVL0
          Write(6,"(a,F15.8)") "Beta_(j+1)",Beta_jp1
          Write(6,"(a,F15.8)") "Beta_(j+1)",Gama_jp1
#endif

C v_(j) --> v_(j-1) and w_(j) --> w_(j-1)

          Call Dcopy(Nsize,Vecs(1,3),1,Vecs(1,1),1)
          Call Dcopy(Nsize,Vecs(1,4),1,Vecs(1,2),1)
          
C v_(j+1) --> v_(j) and  w_(j+1) --> w_(j)

          Call Dcopy(Nsize,Vecs(1,5),1,Vecs(1,3),1)
          Call Dcopy(Nsize,Vecs(1,6),1,Vecs(1,4),1)

          J = J + 1
          Retlabel = 2

C Form the tridiagonal matrix's elements (iteration no. 1)

          Trd(1,1) = -Alpa_jrl 
CSSS          Trd(2,1) = Beta_jp1
CSSS          Trd(1,2) = Gama_jp1

          Goto 20

      Else 

          Alpa_jrl = -Ddot(Nsize,Vecs(1,5),1,Vecs(1,4),1)

#ifdef _DEBUG_LVLM          
          Write(6,"(a,F15.8)") "<L_old|R_new> = ", -Alpa_jrl
#endif

C p_(j) = V_(j+1) - alpa_(j)V_(j)-beta_(j)v_(j-1)

#ifdef _DEBUG_LVL0
          Write(6,"(a,F15.8)") "Beta_(j-1)",Beta_jrl
          Write(6,"(a,F15.8)") "Gama_(j-1)",Gama_jrl
#endif
C Form the tridiagonal matrix's elements (iteration no. j)

          Trd(J,J)   = -Alpa_jrl
          Trd(J-1,J) = -Beta_jrl
          Trd(J,j-1) = -Gama_jrl 

          Call Daxpy(Nsize,Alpa_jrl,Vecs(1,3),1,Vecs(1,5),1)
          Call Daxpy(Nsize,Gama_jrl,Vecs(1,1),1,Vecs(1,5),1)

C s_(j) = w_(j+1) - alpa_(j)w_(j)-beta_(j)w_(j-1)

          Call Daxpy(Nsize,Alpa_jrl,Vecs(1,4),1,Vecs(1,6),1)
          Call Daxpy(Nsize,Beta_jrl,Vecs(1,2),1,Vecs(1,6),1)

          Ps_mod   = Ddot(Nsize,Vecs(1,5),1,Vecs(1,6),1)

          Beta_jp1 = Dsqrt(Dabs(Ps_mod))
          Gama_jp1 = (Ps_mod/Dabs(Ps_mod)) * Beta_jp1

          If (Dabs(Beta_jp1) .Lt. Threshold) Then

            Write(6,"(a,a,i5)") "The <L_(n+1)|R_(n+1)> j=1,...n", 
     +      " is zero at n = ",j
            Write(6,"(a)") "Look ahead algorithm is not implemented"
            Call Errex
          Endif 

          Beta_jp1_inv = 1.0D0/Beta_jp1
          Gama_jp1_inv = 1.0D0/Gama_jp1
 
          Call Dscal(Nsize,Beta_jp1_inv,Vecs(1,5),1)
          Call Dscal(Nsize,Gama_jp1_inv,Vecs(1,6),1)
C
          Beta_jrl = -Beta_jp1
          Gama_jrl = -Gama_jp1

#ifdef _DEBUG_LVL0
          Write(6,"(a,F15.8)") "Beta_(j+1)",Beta_jp1
          Write(6,"(a,F15.8)") "Gama_(j+1)",Gama_jp1
#endif

C v_(j) --> v_(j-1), w_(j) --> w_(j-1)

          Call Dcopy(Nsize,Vecs(1,3),1,Vecs(1,1),1)
          Call Dcopy(Nsize,Vecs(1,4),1,Vecs(1,2),1)

C v_(j+1) --> v_(j), w_(j+1) --> w_(j)

          Call Dcopy(Nsize,Vecs(1,5),1,Vecs(1,3),1)
          Call Dcopy(Nsize,Vecs(1,6),1,Vecs(1,4),1)

C Form the tridiagonal matrix's elements (iteration no. j)
CSSS          Trd(J,J)   = -Alpa_jrl
CSSS          Trd(J+1,J) =  Beta_jp1
CSSS          Trd(J,j+1) =  Gama_jp1

          J = J + 1
          Retlabel = 2 
          Goto 20

       Endif 
      
      Return 
      End

    
