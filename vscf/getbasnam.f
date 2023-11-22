      SUBROUTINE GETBASNAM(IATOM,IZ,NAME)

      IMPLICIT NONE
      INTEGER IATOM,IZ
      CHARACTER*80 NAME
      CHARACTER*80 WRK
      LOGICAL YESNO,SPECIAL
      INTEGER ILINES,IBUB,INDEX,LUZ,I,LINENUM,LENLINE,FNBLNK
      INTEGER IPOS,ICOMMA,IBLNK,LENBASNAM
C
      do 2 i=1,80
         name(i:i) = ' '
    2 continue
C
      LUZ=71
C
C Open the ZMAT file.
C
      INQUIRE(FILE='ZMAT',OPENED=YESNO)
      IF(.NOT.YESNO)THEN
        OPEN(UNIT=LUZ,FILE='ZMAT',FORM='FORMATTED',STATUS='OLD')
      ENDIF
      REWIND(LUZ)
C
C SKIP TITLE LINE AS WELL AS ALL LINES CONCERNING FILE-HANDLING
C
  102 READ(LUZ,'(A)',END=5400)WRK
      IF(WRK(1:1).EQ.'%'.OR.WRK(1:1).EQ.' ') GO TO 102
      WRK = ' '
    1 CONTINUE
      IF(WRK(1:6).NE.'*CRAPS'.AND.WRK(1:6).NE.'*ACES2') then
        READ(LUZ,'(A)',END=5400)WRK
        GOTO 1
      ENDIF
      BACKSPACE(LUZ)
C
C Find out the number of lines of keyword input.
C At this point we are about to read the first line of the ACES2
C namelist.
C
      ILINES=0
   15 CONTINUE
      READ(LUZ,'(A)',END=8012)WRK
      IF(WRK(1:1) .NE. ' ')THEN
        ILINES=ILINES+1
c       READ(LUZ,'(A)',END=8012)WRK
        GOTO 15
      ENDIF
 8012 CONTINUE
C
C At this point we have read the line after the namelist. In other
C words, the current line is the SECOND line after the namelist.
C
      BACKSPACE(LUZ)
C
C The current line is the first after the namelist.
C NOW RETURN TO FIRST LINE OF KEYWORDS.
C
      DO 100 I=1,ILINES
         BACKSPACE(LUZ)
  100 CONTINUE
C
C The current line is the first line of the namelist.
C
      do 200 i=1,ilines
         read(LUZ,'(A)') WRK
  200 continue 
C
      DO 300 I=1,ILINES
         BACKSPACE(LUZ)
  300 CONTINUE
C
      LINENUM=0
  301 CONTINUE
      LINENUM=LINENUM+1
      READ(LUZ,'(A)') WRK
      IBUB=INDEX(WRK,'BASIS=')
      IF(IBUB.EQ.0 .AND. LINENUM.GT.ILINES)THEN
       write(6,*) ' BASIS= string not found '
       call errex
      ENDIF
      IF(IBUB.EQ.0) GOTO 301
C
C Find the starting position of the BASIS= string.
C
      IPOS=INDEX(WRK,'BASIS=')
C
C     Do we have BASIS=SPECIAL?
C
      IBUB = INDEX(WRK,'BASIS=SPECIAL')
      IF(IBUB .NE. 0)THEN
         SPECIAL = .TRUE.
      ENDIF
C
C-----------------------------------------------------------------------
C     BASIS=SPECIAL
C-----------------------------------------------------------------------
      IF(SPECIAL)THEN
C
  302  CONTINUE
       READ(LUZ,'(A)') WRK
       IBUB = INDEX(WRK,' ')
       IF(WRK(1:1) .NE. ' ') GOTO 302
       DO 303 I = 1,IATOM
       READ(LUZ,'(A)') WRK
  303  CONTINUE
C
       IBLNK=fnblnk(wrk)
C
       IF(WRK(2:2) .NE. ':' .AND. WRK(3:3) .NE. ':')THEN
        write(6,*)
     &   ' @GETBASNAM-F, Special basis set. Something very wrong.'
        call errex
       ENDIF
C
       IF(WRK(2:2) .EQ. ':')THEN
        lenbasnam=iblnk-1-1
        do 304 i=1,lenbasnam
        name(i:i) = wrk(2+i:2+i)
  304   continue
       ENDIF
C
       IF(WRK(3:3) .EQ. ':')THEN
        lenbasnam=iblnk-2-1
        do 305 i=1,lenbasnam
        name(i:i) = wrk(3+i:3+i)
  305   continue
       ENDIF
C
       IF(.NOT.YESNO) CLOSE(UNIT=LUZ,STATUS='KEEP')
       return
C
      ENDIF
C
C-----------------------------------------------------------------------
C *** NOT *** BASIS=SPECIAL, Read the basis name from the JOBARC 
C-----------------------------------------------------------------------
C
      CALL GETREC (20, 'JOBARC', 'BASNAMLN', 1, LENBASNAM)
      CALL GETCREC(20, 'JOBARC', 'BASISNAM', LENBASNAM,
     &             NAME)
C
      IF(.NOT.YESNO) CLOSE(UNIT=LUZ,STATUS='KEEP')
C
      RETURN
 5400 WRITE(6,5401)
 5401 FORMAT(T3,'@GTFLGS-F, *ACES2 namelist not found on ZMAT.')
      CLOSE(UNIT=LUZ,STATUS='KEEP')
      CALL ERREX
      END
