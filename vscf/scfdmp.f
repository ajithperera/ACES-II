      SUBROUTINE SCFDMP(FNEW,FOLD,DAMP,LEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FNEW(LEN),FOLD(LEN)
C
      DO 10 I=1,LEN
      FNEW(I) = (FNEW(I) + DAMP * FOLD(I)) / (1.0D+00 + DAMP)
   10 CONTINUE
      RETURN
      END