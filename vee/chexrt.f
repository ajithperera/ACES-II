      SUBROUTINE CHEXRT(NROOTS)
C
C SG 2/11/97
C This subroutine looks for the string
C     '%follow root*'
C If it sees this string, it will read from the next line a single
C  integer.  That will be the number of roots that will calculated.
C  The primary use for this subroutine is in doing P-EOM-MBPT(2)
C  gradient calculations, where the user can specify a large number
C  for ESTATE_SYM, but then after the block Davidson roots have been
C  reordered, he can follow whichever root he wants.
C
C      IMPLICIT NONE
      INTEGER NROOTS,NREAD
C
      INTEGER INDEX
      CHARACTER*80 STRING
C
      OPEN(UNIT=30,FILE='ZMAT',FORM='FORMATTED')
      REWIND(30)
C
 1    CONTINUE
      READ(30,'(A)',END=800)STRING
      IF (INDEX(STRING,'%FOLLOW ROOT*').EQ.0 .AND.
     &   INDEX(STRING,'%follow root*').EQ.0) THEN
        GOTO 1
      ELSE
        READ(30,*,END=900,ERR=950)NREAD
        IF (NREAD .LE. 0) THEN
          WRITE(6,1000)
 1000     FORMAT(T3,'@CHEXRT-F, Input for root to be followed must',
     &       ' be positive.')
          CLOSE(UNIT=30,STATUS='KEEP')
          CALL ERREX
        ELSEIF (NREAD .GT. NROOTS) THEN
          WRITE(6,1010) NREAD,NROOTS
 1010     FORMAT(T3,'@CHEXRT-F, Root followed',I5,' must be less than',
     &       ' total number of roots',I5,'.')
          CLOSE(UNIT=30,STATUS='KEEP')
          CALL ERREX
        ELSE
          NROOTS = NREAD
        ENDIF
      ENDIF
 800  CONTINUE
      CLOSE(UNIT=30,STATUS='KEEP')
      RETURN
C
 900  CONTINUE
      WRITE(6,2000)
 2000 FORMAT(T3,'@CHEXRT-F, Root to follow not specified.')
      CLOSE(UNIT=30,STATUS='KEEP')
      CALL ERREX
C
 950  CONTINUE
      WRITE(6,2010)
 2010 FORMAT(T3,'@CHEXRT-F, Root to follow must be given as a single',
     &   ' integer.')
      CLOSE(UNIT=30,STATUS='KEEP')
      CALL ERREX
C
      RETURN
      END
