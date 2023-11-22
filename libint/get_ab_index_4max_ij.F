      Subroutine Get_ab_index_4max_ij(D_kp1,Nbasis,Nbasis2,Jmax,
     &                                Dmax_old)

      Implicit Double Precision(A-H,O-Z)

      Dimension D_kp1(Nbasis2,Nbasis2)
      Logical Degenerate 

      Data Dzero /0.0D0/
     
      Dmax_new   = Dzero
      Dmax_old   = Dzero
      Degenerate = .False.
 
      Do J = 1, Nbasis2
         Do I = 1, Nbasis2
            Dmax_new = Max(Dmax_new,D_kp1(I,J))
            If (Dmax_new .Gt. Dmax_old) Then
                Imax = I
                Jmax = J
                Dmax_old = Dmax_new
                Dmax_new = Dzero
            Elseif (Dmax_new .Eq. Dmax_old) Then 
                Degenerate = .True.
            Endif  
         Enddo
      Enddo 

      If (Degenerate) Then
         Write(6,"(2a)") " One or more degenerte maximum values are",
     &                   " present in the residual."
         Write(6,"(a)") " First found first used criteria is used."
      Endif 

      Return
      End
            
        
