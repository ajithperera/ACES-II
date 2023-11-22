      SUBROUTINE RELVEC(R,E,C1,C2)
C-----------------------------------------------------------------------
C     Routine to calculate internuclear distance and unit vectors for
C     internuclear vector in molecular frame.
C     Modernized version of routine given by Pople and Beveridge.
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION R,E,C1,C2
      DOUBLE PRECISION X
      INTEGER I
      DIMENSION E(3),C1(3),C2(3)
C
      X=0.0D+00
      DO 10 I=1,3
      E(I) = C2(I)-C1(I)
      X    = X + E(I)**2
   10 CONTINUE
C
      R = DSQRT(X)
C
      IF(R.GT.0.000001D+00)THEN
       DO 20 I=1,3
        E(I) = E(I)/R
   20  CONTINUE
      ENDIF
C
      RETURN
      END
