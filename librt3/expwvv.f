      SUBROUTINE EXPWVV(T2,W1,W2,NAB,NA,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,AB
      DIMENSION T2(NAB),W1(NA,NB),W2(NB,NA)
      DO   20 B=1,NB
      DO   10 A=1,NA
      AB = (B-1) * NA + A
      T2(AB) = T2(AB) + W1(A,B) - W2(B,A)
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
