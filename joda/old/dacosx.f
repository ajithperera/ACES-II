      DOUBLE PRECISION FUNCTION DACOSX(VALUE)
C
C FUNCTION WHICH ROUNDS ARGUMENTS VERY NEAR -1.0/1.0 TO -1.0/1.0.  
C NECESSARY FOR SOME THINGS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (TOL = 1.D-10)
      DATA ONE /1.0/
      X=ABS(ABS(VALUE)-ONE)
      IF(X.LT.TOL)VALUE=SIGN(ONE,VALUE)
      DACOSX=DACOS(VALUE)
      RETURN
      END
