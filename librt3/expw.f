      SUBROUTINE EXPW(T3,W1,W2,NAB,NC,NA,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,AB
      DIMENSION T3(NAB,NC),W1(NA,NB,NC),W2(NB,NA,NC)
      DO   30 C=1,NC
      DO   20 B=1,NB
      DO   10 A=1,NA
      AB = (B-1) * NA + A
      T3(AB,C) = T3(AB,C) + W1(A,B,C) - W2(B,A,C)
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
      RETURN
      END