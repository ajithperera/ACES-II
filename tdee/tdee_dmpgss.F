      SUBROUTINE TDEE_DMPGSS(R, RLEN, NAME)
C
C This routine dumps the R vector into a file to be picked up and used
C as a guess in a later CC calculation.
C
C R contains the R amplitudes in the following order
C R(AA-BB-AAAA-BBBB-ABABA) for UHF 
C R(AA-BB-ABABA) for RHF 
C
      IMPLICIT NONE
C
      INTEGER RLEN
      DOUBLE PRECISION R(RLEN)
      CHARACTER *8 NAME
C
      INTEGER NAMLEN
      LOGICAL YESNO
      CHARACTER *80 FULNAM
C
      CALL GFNAME(NAME, FULNAM, NAMLEN)
C
      INQUIRE(FILE=FULNAM(1:NAMLEN), EXIST=YESNO)
      IF (YESNO) THEN
        OPEN(UNIT=94, FILE=FULNAM(1:NAMLEN), STATUS='OLD',
     &     FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
        CLOSE(UNIT=94, STATUS='DELETE')
      ENDIF
C
      OPEN(UNIT=94, FILE=FULNAM(1:NAMLEN), STATUS='NEW',
     &   FORM='UNFORMATTED', ACCESS='SEQUENTIAL')
      WRITE(94) R
      CLOSE(UNIT=94, STATUS='KEEP')
C
      RETURN
      END
