            SUBROUTINE MATPRT2(A,NR,NC,LM2,LM3,C)
C     *
C     PRINT THE MATRIX A(LM2,LM3), NR ROWS AND NC COLUMNS.
C     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(LM2,LM3),C(NR,NR)
      KA     = 1
      KC     = 10
   10 KB     = MIN(KC,NC)
      WRITE(6,500) (I,I=KA,KB)
      N      = 0
      KB=KB+1
      WRITE(6,520) (C(I,I),I=KB,KA,-1)
      PRINT*,''
      DO 20 I=1,NR
      WRITE(6,520) I,(A(I,J),J=KA,KB)
      N      = N+1
      IF(N.LT.10) GO TO 20
      WRITE(6,510)
      N      = 0
   20 CONTINUE
      KB=KB-1
      IF(KB.EQ.NC) RETURN
      KA     = KC+1
      KC     = KC+10
      GO TO 10
  500 FORMAT(// 9X,I5,9I12)
  510 FORMAT(   1X)
  520 FORMAT(   1X,I4,10F12.5)
      END
