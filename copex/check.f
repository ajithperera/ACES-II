      SUBROUTINE CHECK(A,N,B)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(N)
C
      B=0.D+0
C
      DO 1 I=1, N
         B=B+A(I)*A(I)
 1    CONTINUE
C
      RETURN
      END
