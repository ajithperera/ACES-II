      Subroutine Massw_cfc(Cfc_cart,Atmmass,Ndim,Nmodes,Nreals)

      Implicit Double Precision(A-H,O-Z)

      Dimension Cfc_cart(Ndim,Ndim,Nmodes)
      Dimension Atmmass(Nreals)

      Data Dnull, Thres /0.0D0,1.0D-09/

      Do I = 1, Nmodes
         Do J = 1, Ndim
            Do K = 1, Ndim
               Fact = Dsqrt(Atmmass(1+(J-1)/3)*Atmmass(1+(K-1)/3))
               If (Fact .Lt. Thres) Then
                  Cfc_cart(K,J,I) = Dnull
               Else
                  Cfc_cart(K,J,I) = Cfc_cart(K,J,I)/Fact
               Endif
            Enddo
         Enddo
      Enddo
    
      Return 
      End
