      LOGICAL FUNCTION QSAME(Q1,Q2)
C
C THIS ROUTINE DETERMINES IF THE COORDINATE VECTORS Q1 AND Q2 ARE
C  THE SAME
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Q1(3),Q2(3)
      QSAME=.TRUE.
      DO 10 I=1,3
       Z=ABS(Q1(I)-Q2(I))
       IF(Z.GT.1.D-8)QSAME=.FALSE.
10    CONTINUE
      RETURN
      END
