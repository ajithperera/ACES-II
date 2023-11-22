      SUBROUTINE AB2SD(DENSA,DENSB,DENSSD,LEN)
C
C  This routine constructs the singles and doubles density matrices from the
C  alpha and beta density matrices.  This is done using the formula:
C
C
C      D(doubles) = 2*D(beta)
C
C      D(singles) = D(alpha) - D(beta)
C
CEND
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION DENSA(LEN),DENSB(LEN),DENSSD(2*LEN)
      COMMON /POPUL/ NOCCA(8),NOCCB(8)
#include "symm2.com"
C
      DATA ONEM /-1.0/
      DATA TWO /2.0/
C
      DO 100 I=1,NIRREP
        IF(NBFIRR(I).EQ.0) GOTO 100
        IF(NOCCA(I).EQ.NOCCB(I)) THEN
          CALL SCOPY(ITRILN(I),DENSB(ITRIOF(I)),1,DENSSD(ITRIOF(I)),1)
          CALL SSCAL(ITRILN(I),TWO,DENSSD(ITRIOF(I)),1)
          CALL ZERO(DENSSD(LEN+ITRIOF(I)),ITRILN(I))
        ELSE
          CALL SCOPY(ITRILN(I),DENSA(ITRIOF(I)),1,
     &               DENSSD(LEN+ITRIOF(I)),1)
          CALL SAXPY(ITRILN(I),ONEM,DENSB(ITRIOF(I)),1,
     &               DENSSD(LEN+ITRIOF(I)),1)
        ENDIF
  100 CONTINUE
C
      RETURN
      END
