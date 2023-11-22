      Subroutine Lanczos_light_v0(Nsize,Vecs,Info)

      Implicit Double Precision (A-H,O-Z)
      Integer Retlabel, Revcom, Info(4)

      Parameter (Threshold = 1.0D-10)

      Dimension Vecs(Nsize,6)

      Common /Beta/ Beta_jrl, J
      Data Retlabel /1/

      Revcom  = Info(2)
      Info(2) = 0
      If (Revcom .EQ. 0) J = 1

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "--- Entered Lanczos_light ---"
      Write(6,"(a,i4,1x,i2)") "The counter,retlabel :", j,retlabel
#endif

      Goto (10,30,40) Retlabel 
      
 10   Continue 

      Call Tdee_normaliz_init_vec(Vecs(1,1),Nsize)
      RL_Norm = Ddot(Nsize,Vecs(1,1),1,Vecs(1,2),1) 
      Call Dscal(Nsize,1.0D0/RL_NORM,Vecs(1,2),1)

#ifdef _DEBUG_LVL0
      RL_Norm = Ddot(Nsize,Vecs(1,1),1,Vecs(1,2),1) 
      Write(6,"(a,F15.8)")  "@-Entry; RL_Norm: ", RL_Norm
#endif 

C Check whether <R|L> = 1.0 for the starting vectors
     
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
      Write(6,"(a)") "Returning 1"
      Return

 30   Continue 

C Do <w_j+1| = <w_j| A

      Info(2)  = 2
      Info(3)  = 4
      Info(4)  = 6
      Retlabel = 3
      Write(6,"(a)") "Returning 2"
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

          Beta_jp1 = Dsqrt(Ddot(Nsize,Vecs(1,5),1,Vecs(1,6),1))
          beta_jp1_inv = 1.0D0/Beta_jp1
      
          Call Dscal(Nsize,beta_jp1_inv,Vecs(1,5),1)
          Call Dscal(Nsize,beta_jp1_inv,Vecs(1,6),1)

          Beta_jrl = -Beta_jp1
#ifdef _DEBUG_LVL0
          Write(6,"(a,F15.8)") "Beta_(j+1)",Beta_jp1
#endif

C v_(j) --> v_(j-1) and w_(j) --> w_(j-1)

          Call Dcopy(Nsize,Vecs(1,3),1,Vecs(1,1),1)
          Call Dcopy(Nsize,Vecs(1,4),1,Vecs(1,2),1)
          
C v_(j+1) --> v_(j) and  w_(j+1) --> w_(j)

          Call Dcopy(Nsize,Vecs(1,5),1,Vecs(1,3),1)
          Call Dcopy(Nsize,Vecs(1,6),1,Vecs(1,4),1)

          J = J + 1
          Retlabel = 2
          Goto 20

      Else 

          Alpa_jrl = -Ddot(Nsize,Vecs(1,5),1,Vecs(1,4),1)
          
          Write(6,"(a,F15.8)") "<L_old|R_new> = ", -Alpa_jrl

C p_(j) = V_(j+1) - alpa_(j)V_(j)-beta_(j)v_(j-1)

#ifdef _DEBUG_LVL0
           
          Write(6,"(a,F15.8)") "Beta_(j-1)",Beta_jrl
#endif
          Call Daxpy(Nsize,Alpa_jrl,Vecs(1,3),1,Vecs(1,5),1)
          Call Daxpy(Nsize,Beta_jrl,Vecs(1,1),1,Vecs(1,5),1)

C s_(j) = w_(j+1) - alpa_(j)w_(j)-beta_(j)w_(j-1)

          Call Daxpy(Nsize,Alpa_jrl,Vecs(1,4),1,Vecs(1,6),1)
          Call Daxpy(Nsize,Beta_jrl,Vecs(1,2),1,Vecs(1,6),1)

          Beta_jp1 = Dsqrt(Ddot(Nsize,Vecs(1,5),1,Vecs(1,6),1))

          If (Dabs(Beta_jp1) .Lt. Threshold) Then

            Write(6,"(a,a,i5)") "The <L_(n+1)|R_(n+1)> j=1,...n", 
     +      " is zero at n = ",j
            Write(6,"(a)") "Look ahead algorithm is not implemented"
            Call Errex
          Endif 

          beta_jp1_inv = 1.0D0/Beta_jp1

          Call Dscal(Nsize,beta_jp1_inv,Vecs(1,5),1)
          Call Dscal(Nsize,beta_jp1_inv,Vecs(1,6),1)
C
          Beta_jrl = -Beta_jp1

#ifdef _DEBUG_LVL0
          Write(6,"(a,F15.8)") "Beta_(j+1)",Beta_jp1
#endif

C v_(j) --> v_(j-1), w_(j) --> w_(j-1)

          Call Dcopy(Nsize,Vecs(1,3),1,Vecs(1,1),1)
          Call Dcopy(Nsize,Vecs(1,4),1,Vecs(1,2),1)

C v_(j+1) --> v_(j), w_(j+1) --> w_(j)

          Call Dcopy(Nsize,Vecs(1,5),1,Vecs(1,3),1)
          Call Dcopy(Nsize,Vecs(1,6),1,Vecs(1,4),1)

          J = J + 1
          Retlabel = 2 
          Goto 20

       Endif 
      
      Return 
      End

    
