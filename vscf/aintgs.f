      SUBROUTINE AINTGS(A,RHO,KMAX)
C-----------------------------------------------------------------------
C     Routine for evaluating auxiliary A function up to K.
C     This code is based on the code given by Pople and Beveridge.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION A,RHO
      INTEGER KMAX,K
      DIMENSION A(0:16)
C
      A(0) = DEXP(-RHO)/RHO
C
      DO 10 K=1,KMAX
      A(K) = (A(K-1)*DFLOAT(K) + DEXP(-RHO))/RHO
   10 CONTINUE
      RETURN
      END
