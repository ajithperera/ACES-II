      SUBROUTINE DTABLE(D)
      IMPLICIT NONE
      DOUBLE PRECISION D,A,B,FACT
      INTEGER LA,LB,M
      DIMENSION D(0:3,0:3,0:3)
C
      CALL ZERO(D,4*4*4)
C
      DO 50 LB=0,3
      DO 40 LA=0,3
      DO 30  M=0,MIN(LA,LB)
C
      A = 0.5D+00 * DFLOAT(2*LA+1) * FACT(LA-M) / FACT(LA+M)
      B = 0.5D+00 * DFLOAT(2*LB+1) * FACT(LB-M) / FACT(LB+M)
C
      D(M,LA,LB) = (FACT(M+1)/8.0D+00)**2 * DSQRT(A) * DSQRT(B)
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
      RETURN
      END
