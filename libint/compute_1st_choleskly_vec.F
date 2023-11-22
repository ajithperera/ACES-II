      Subroutine Compute_1st_cholesky_vec(Ints,C_k1,Imax,Jmax,Aa1,Aa2,
     &                                    Bb1,Bb2,Ijij_norm,Nbasis,
     &                                    Ncfps_a,Ncfps_b,Den_tol)

      Implicit Double Precision(A-H,O-Z)
    
      Double Precision Ints, Ijij_norm
      Integer Aa1,Aa2,Bb1,Bb2

      Dimension Ints(1:Nbasis,1:Nbasis,Aa1:Aa2,Bb1:Bb2)
      Dimension C_k1(1:Nbasis,1:Nbasis,1:1) 

      If (Ijij_norm .Le. Den_tol) Then
         Write(6,"(a,ES7.2E2,a)") "The norm is smaller than ",Tol,
     &               " and the division is numerically unstable."
         Call Errex
      Endif 
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The <mu nu|1/r12|lam 1> integral matrix",IJIJ_norm
      do a = aa1,aa2
      do b = bb1,bb2
      Write(6, "(5(1x,F15.10))") ((Ints(d,c,b,a),c=1,nbasis),
     &                             d=1,nbasis)
      enddo
      enddo
#endif 

      Do J = 1, Nbasis 
         Do I = 1, Nbasis 
            C_k1(I,J,1) = Ints(I,J,Imax,Jmax)/Ijij_norm 
         Enddo
      Enddo 

#ifdef _DEBUG_LVL0
      Print*, "1st Cholesky vector"
      call output(C_k1,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 
            
      Return 
      End


