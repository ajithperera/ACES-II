      SUBROUTINE TDEE_REDGSS(T, TLEN, NAME)
C
C This routine reads in a T vector written out by DMPTGSS.  This T vector can
C be used as an initial guess to a CC calculation.
C R contains the R amplitudes in the following order
C R(AA-BB-AAAA-BBBB-ABABA) for UHF
C R(AA-BB-ABABA) for RHF
C
      IMPLICIT NONE
C
      INTEGER TLEN
      DOUBLE PRECISION T(TLEN)
      LOGICAL EXISTS
      CHARACTER *8 NAME
C
      INTEGER NAMLEN
      CHARACTER *80 FULNAM
C
      CALL GFNAME(NAME, FULNAM, NAMLEN)
      INQUIRE(FILE=FULNAM(1:NAMLEN), EXIST=EXISTS)
C
      IF (EXISTS) THEN
        OPEN(UNIT=94, FILE=FULNAM(1:NAMLEN), STATUS='OLD',
     &     FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
        READ(94) T
        CLOSE(UNIT=94, STATUS='KEEP')
      ENDIF
C
      RETURN
      END
