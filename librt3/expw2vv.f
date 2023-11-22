      SUBROUTINE EXPW2VV(T2,W1,NAB,NA,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,AB
      DIMENSION T2(NAB),W1(NA,NB)
      INDEX(I) = I * (I-1) / 2
      DO   20 B=2,NB
      DO   10 A=1,B-1
      AB = INDEX(B-1) + A
      T2(AB) = T2(AB) + W1(A,B) - W1(B,A)
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
