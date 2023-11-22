      SUBROUTINE RDTGSS(T, TLEN, NAME, EXISTS)
C
C This routine reads in a T vector written out by DMPTGSS.  This T vector can
C be used as an initial guess to a CC calculation.
C
C T contains the T amplitudes in the following order
C    T1AI   (List 90,1)
C    T1ai   (List 90,2) (if IUHF .NE. 0)
C    T2abij (List 45)   (if IUHF .NE. 0)
C    T2ABIJ (List 44)
C    T2AbIj (List 46)
C
C SG 7/22/98
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
