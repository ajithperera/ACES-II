      Subroutine Vecprd_22(Nrow_a,Nrow_b,Ncol_i,Ncol_j,Lvo,Dvo,
     +                     Target,Fact)

      Implicit Double Precision (A-H, O-Z)

      Dimension Target(Nrow_a,Ncol_i,Nrow_b,Ncol_j),
     +          Dvo(Nrow_b,Ncol_j)
      Double Precision Lvo(Nrow_a,Ncol_i)

      Do L_c = 1, Ncol_i
      Do L_r = 1, Nrow_a
  
         Do M_c = 1, Ncol_j
         Do M_r = 1, Nrow_b

            Target(L_r,L_c,M_r,M_c) = Target(L_r,L_c,M_r,M_c) + 
     +                                Lvo(L_r,L_c)*Dvo(M_r,M_c)*Fact 
         Enddo
         Enddo 

      Enddo
      Enddo 

      Return
      End 
