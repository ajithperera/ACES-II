      Subroutine Update_diag_residuals(D_km1,D_k,C_k,Mu_max,Nu_max,
     &                                 Nbasis,Aa1,Aa2,Bb1,Bb2,Kmax,
     &                                 Lmax,Ijij_norm,Dmax_old,K)

      Implicit Double Precision(A-H,O-Z)
      Integer Aa1,Aa2,Bb1,Bb2

      Dimension D_km1(Nbasis,Nbasis)
      Dimension D_k(Nbasis,Nbasis)
      Dimension C_k(Nbasis,Nbasis,K)
      Dimension Mu_max(Nbasis)
      Dimension Nu_max(Nbasis)
      Double Precision Ijij_norm 
      Logical Degenerate

      Data Dzero /0.0D0/

      Dmax_old  = Dzero
      Ijij_norm = Dzero 

#ifdef _DEBUG_LVL0
      Print*, "@-update_diag D_km1",K
      call output(D_km1,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      Print*, "@-update_diag ck 1->K" 
      do i=1,k
      call output(C_k(1,1,i),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
      enddo 
#endif 

      Do I = 1, Nbasis 
         Do J = 1, Nbasis 
            D_k(J,I) = D_km1(J,I) - C_k(J,I,K)*C_k(J,I,K)
            D = D_k(J,I) 
            Print*, D
            Ijij_norm = Ijij_norm + D*D
            Dmax_new = D
            If (Dmax_new .Gt. Dmax_old) Then
               Kmax = I
               Lmax = J
               Dmax_old = Dmax_new
            Elseif (Dmax_new .Eq. Dmax_old) Then
               Degenerate = .True.
            Endif
         Enddo 
      Enddo 

      Ijij_norm = Dsqrt(Dmax_old)
      Mu_max(K) = Lmax
      Nu_max(K) = Kmax

      Call Dcopy(Nbasis*Nbasis,D_k,1,D_km1,1)

#ifdef _DEBUG_LVL0
      Print*
      Print*, "Kmax and Lmax", kmax,Lmax
      Print*, "Dmax", Dmax_old
      Print*
      Write(6,"(6(1x,i3))") (Mu_max(i),i=1,K)
      Write(6,"(6(1x,i3))") (Nu_max(i),i=1,K)
      Print*, "@-update_diag D_k",K
      call output(D_k,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

      If (Degenerate) Then
         Write(6,"(2a)") " One or more degenerte maximum values are",
     &                   " present in the residual."
         Write(6,"(a)") " First found first used criteria is used."
      Endif

      Return
      end 
