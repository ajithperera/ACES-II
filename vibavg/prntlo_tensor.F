      SUBROUTINE PRNTLO_TENSOR(NOUT, A, B, N, M, NNN, MMM)
C
C This routine prints a lower-trangle of a square-matrix in 
C neat format. The parameter WID controls the number of columns
C printed
C
      IMPLICIT REAL*8(A-H, O-Z)
C
      DIMENSION A(NNN,MMM)
      CHARACTER*6 B(NNN)
C
      INTEGER WID
C     
      PARAMETER(WID = 5)
 1    FORMAT(2X, A6, 1X, 5F12.4)
 12   FORMAT(6X, 5(8X,A6))
C
 11   FORMAT(/)
      MM = M/WID
      JCOUNT = 1
C
      IF (MM .EQ. 0) GO TO 6
C
      DO 10 II = 1, MM 
C
         JP = (II - 1)*WID + 1
         JK = II*WID
C
         WRITE(NOUT, 11)
         WRITE(NOUT, 12) (B(I), I = JP, JK)
         WRITE(NOUT, 11)
C
         DO 20 I = 1, N
C
            WRITE(NOUT, 1) B(I), (A(I,J), J = JCOUNT, MIN(I - 1, 
     &                           JCOUNT + 4))
C
 20      CONTINUE
C
         JCOUNT = JCOUNT + 5
C
 10   CONTINUE
C
 6    CONTINUE
C
      MA = MM*WID + 1
      IF (MA .GT. M) RETURN
C
      WRITE(NOUT, 11)
      WRITE(NOUT, 12)(B(I), I = MA, M)
      WRITE(NOUT, 11)
C
      DO 5 I = 1, N
         WRITE(NOUT, 1) B(I),(A(I,J),J = MA, I-1)
    5 CONTINUE
C
      RETURN
      END
