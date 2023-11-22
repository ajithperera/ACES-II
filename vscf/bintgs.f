      SUBROUTINE BINTGS(B,RHO,KMAX)
C-----------------------------------------------------------------------
C     Routine for evaluating auxiliary B function up to KMAX.
C     This code is based on the code given by Pople and Beveridge.
C     Currently we are hoping we can always use the exponential
C     formula (except for 0).
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION B,RHO,TOL,EXPRHO,EXPMRHO
      INTEGER KMAX,K
      DATA TOL /0.000001D+00/
      DIMENSION B(0:16)
C
      IF(DABS(RHO) .LT. TOL)THEN
       DO 10 K=0,KMAX
       B(K) = (1.0D+00 + (-1.0D+00)**K) / DFLOAT(K+1)
   10  CONTINUE
      ELSE
C
       EXPRHO  = DEXP(RHO)
       EXPMRHO = 1.0D+00/EXPRHO
C
       B(0) = (EXPRHO - EXPMRHO)/RHO
       DO 20 K=1,KMAX
       B(K) = ( B(K-1)*DFLOAT(K) + (-1.0D+00)**K*EXPRHO - EXPMRHO )/RHO
   20  CONTINUE
C
      ENDIF
C
      write(6,*) ' bintgs ',rho,kmax
      write(6,*) b
C
      RETURN
      END
