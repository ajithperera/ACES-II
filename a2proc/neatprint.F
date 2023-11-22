C
      SUBROUTINE NEATPRINT(NOUT, A, N, M, NNN, MMM)
C
C This routine prints a sqaure matrix in a neat format.
C The parameter WID controls the number of columns 
C printed. 
C
C  NOUT     : The unit number to which the output being 
C             written.
C  A        : Input matrix.
C  N,M      : Actual dimensions of the matrix.
C  NNN,MMM  : Logical dimensions
C  
      IMPLICIT REAL*8(A-H,O-Z)
C     
      DIMENSION A(NNN, MMM)
      INTEGER WID
C
      PARAMETER(WID = 5)
 110  FORMAT(1X, I3, 2X, 5F15.9)
 120  FORMAT(4X, 5(8X, I4))
C
 130  FORMAT(/)
      MM = M/WID
      IF (MM .EQ. 0) GO TO 100
C
      DO 10 II = 1, MM 
         JP = (II - 1)*WID + 1
         JK = II*WID
C
CSSS         WRITE(NOUT, 130)
         WRITE(NOUT, 120)(I, I = JP, JK)
         WRITE(NOUT, *)
C
         DO 20 I = 1, N
            WRITE(NOUT, 110) I, (A(I,J), J = JP, JK)
 20      CONTINUE
C
 10   CONTINUE
C
 100  CONTINUE
C
      MA = MM*WID + 1
      IF (MA .GT. M) RETURN
C
CSSS      WRITE(NOUT, 130)
      WRITE(NOUT, 120)(I, I = MA, M)
      WRITE(NOUT, *)
      DO 30 I = 1, N
         WRITE(NOUT, 110) I, (A(I,J),J=MA,M)
 30   CONTINUE
C
      RETURN
      END

