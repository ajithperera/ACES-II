      Subroutine Sqrt_diags(H,Nmodes,Max_asym)

      Implicit Double Precision(A-H,O-Z)
      Double Precision Max_asym

      Dimension H(Nmodes,Nmodes)

      Data DNull /0.0D0/

      Max_asym = DNull

      Do J = 1, Nmodes 
         Do I = 1, Nmodes 
            Max_asym = Max(Max_asym,(Dabs(H(I,J))-Dabs(H(J,I))))
         Enddo 
      Enddo 

      Do I = 1, Nmodes
         If (H(I,I) .Le. Dnull) Then
             H(I,I) = Dsqrt(-H(I,I))
         Else
             H(I,I) = Dsqrt(H(I,I))
         Endif 
      Enddo 

      Return
      End 
