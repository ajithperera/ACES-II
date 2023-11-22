      SUBROUTINE RDRGSS(T,TLEN,NAME,EXISTS)
C
C This routine reads in a R vector written out by DMPRGSS.  This R vector can
C be used as an initial guess to a CC calculation.
C
C R contains the R amplitudes in the following order
C    R1AI   (List 490,1)
C    R1ai   (List 490,2) (if IUHF .NE. 0)
C    R2abij (List 445)   (if IUHF .NE. 0)
C    R2ABIJ (List 444)
C    R2AbIj (List 446)
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
      CALL GFNAME(NAME, FULNAM, NAMLEN)
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
