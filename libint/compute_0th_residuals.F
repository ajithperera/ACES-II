      Subroutine Compute_0th_residuals(Ints_current,D_current,Mu_Max,
     +                                 Nu_max,Nbasis,Imax,Jmax,
     +                                 Ijij_norm)

      Implicit Double Precision (A-H, O-Z)

      Double Precision Ints_current 
      Double Precision Ijij_norm
      Double Precision Max_int_new,Max_int_old
      Integer Aa1,Aa2,Bb1,Bb2
      Logical Degenerate 

      Dimension Ints_current(Nbasis,Nbasis,Nbasis,Nbasis)
      Dimension D_current(Nbasis,Nbasis)
      Dimension Mu_max(Nbasis)
      Dimension Nu_max(Nbasis)

      Data Dzero /0.0D0/

      Do I = 1, Nbasis 
         Do J = 1, Nbasis 
            D_current(J,I) = Ints_current(J,I,J,I) 
         Enddo
      End do
#ifdef _DEBUG_LVL0
      write(6,"(a)") "The I(I,J|I,J) ints"
      call output(D_current,1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif 

      Max_int_new = Dzero
      Max_int_old = Dzero
      Ijij_norm   = Dzero
      Do b = 1, Nbasis 
         Do a = 1, Nbasis 
            Ijij_norm = Ijij_norm + D_current(a,b)*
     &                              D_current(a,b)
            Max_int_new = D_current(a,b)
            If (Max_int_new .Gt. Max_int_old) Then
               Imax = a
               Jmax = b
               Max_int_old = Max_int_new
           Elseif (Max_int_new  .Eq. Max_int_old) Then
                Degenerate = .True.
           Endif
         Enddo
      Enddo

      Ijij_norm = Dsqrt(Max_int_old)
      Mu_max(1) = Imax
      Nu_max(1) = Jmax

      If (Degenerate) Then
         Write(6,"(2a)") " One or more degenerte maximum values are",
     &                   " present in the residual."
         Write(6,"(a)") " First found first used criteria is used."
      Endif


      Return
      End

