      subroutine Compute_kp1_cholesky_vec(D_kp1,C_kp1,Nbasis,Nbasis2,
     &                                    Ncfps_a,Ncfps_b,Aa1,Aa2,Bb1,
     &                                    Bb2,Kmax,Lmax,Ijij_norm,
     &                                    Den_tol,K)

      Implicit Double Precision(A-H,O-Z)
    
      Double Precision Ints, Ijij_norm
      Integer Ab,Aa1,Aa2,Bb1,Bb2

      Dimension D_kp1(1:Nbasis2)
      Dimension C_kp1(1:Nbasis2,1:K)

      If (Ijij_norm .Le. Den_tol) Then
         Write(6,"(a,ES7.2E2,a)") "The norm is smaller than ",Tol,
     &               " and the division is numerically unstable."
         Call Errex
      Endif 

      Do J = 1, Nbasis2
         C_kp1(J,K) = D_kp1(J)/Ijij_norm
      Enddo 

#ifdef _DEBUG_LVL0
      Print*, "@_kp1_cholesky_vec: C_kp1",k
      do i=1,k
      Call output(C_kp1(1,i),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      enddo 
#endif 
            
      Return 
      End


