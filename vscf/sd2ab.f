      SUBROUTINE SD2AB(DENSSD,DENSA,DENSB,LEN)
C
C  This subroutine constructs the alpha and beta density matrices from
C  the singles and doubles density matrices in an ROHF calculation. 
C  This is done using the formula:
C
C    D(alpha) = 0.5*D(doubles) + D(singles)
C
C    D(beta)  = 0.5*D(doubles)
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION DENSSD(2*LEN),DENSA(LEN),DENSB(LEN)
C
      DATA HALF /0.5/
      DATA ONE /1.0/
C
C  Copy the doubles density into both the alpha and beta parts and
C  divide by 2.
C
      CALL SCOPY(LEN,DENSSD,1,DENSA,1)
      CALL SCOPY(LEN,DENSSD,1,DENSB,1)
      CALL SSCAL(LEN,HALF,DENSA,1)
      CALL SSCAL(LEN,HALF,DENSB,1)
C
C  Now add the singles density to the alpha part.
C
      CALL SAXPY(LEN,ONE,DENSSD(LEN+1),1,DENSA,1)
C
C  All done!
C
      RETURN
      END
