      FUNCTION FACT(N)
C-----------------------------------------------------------------------
C     Simple factorial function based on Pople and Beveridge.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION FACT,PROD
      INTEGER N,I
      PROD=1.0D+00
      IF(N.GT.0)THEN
       DO 10 I=1,N
       PROD=PROD * DFLOAT(I)
   10  CONTINUE
      ENDIF
      FACT=PROD
      RETURN
      END
