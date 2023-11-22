C
      Subroutine GaussL(X1, X2, X, W, N)
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C &  Added 08/93 by Ajith to calculate diamagnetic contribution  &
C &  to the NMR coupling constant. Calculates Weights and roots  &
C &  for Gauss-Legendre numerical integration. Taken from        &
C &  Numerical Recipes. Given the lower and upper limits of      &
C &  integration X1, X2 and given N this routine returns the     &
C &   weights (W) and abscissas (X) of length N.                 &
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      Implicit Double Precision (A-H, O-Z)
      Dimension X(N), W(N)
      Parameter (EPS = 3.0D-14)
      M = (N +1)/2
      XM = 0.50D00*(X2 + X1)
      XL = 0.50D00*(X2 - X1)
C
      Do 10 I = 1, M
         Z = Cos(3.141592654D0*(I - 0.25D+00)/(N + 0.50D+00))
C 
 100     Continue
         P1 = 1.0D+00
         P2 = 0.0D+00
C
         Do 20 J = 1, N
            P3 = P2
            P2 = P1
            P1 = ((2.0D+00*J - 1.00D+00)*Z*P2-(J - 1.0D+00)*P3)/J
 20      Continue
         PP = N*(Z*P1 - P2)/(Z*Z - 1.0D+00)
         Z1 = Z
         Z = Z1 - P1/PP
         If (Abs(Z - Z1) .gt. EPS) goto 100
         X(I) = XM - XL*Z
         X(N + 1 - I) = XM + XL*Z
         W(I) = 2.0D+00*XL/((1.0D+00 - Z*Z)*PP*PP)
         W(N + 1 - I) = W(I)
 10   Continue
      Return
      End
