      SUBROUTINE FNDMAXD(N, VEC, DIFF, ILOC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VEC(N)
C
      AMAX = 0.D0
      DO 10 I=1, N
         IF (ABS(VEC(I)).GT.AMAX) THEN
            AMAX = ABS(VEC(I))
            ILOC = I
         ENDIF
 10   CONTINUE
C
      DIFF = ABS(VEC(ILOC))
C
      RETURN
      END
