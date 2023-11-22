      SUBROUTINE GETBASINF(IZATOM,BASNAM,BASTYP,NSFULL,NPFULL,NDFULL)
      IMPLICIT NONE
      INTEGER IZATOM,NSFULL,NPFULL,NDFULL
      CHARACTER*80 BASNAM
      CHARACTER*6  BASTYP
      CHARACTER*80 TITLE,WRK
      INTEGER INDEX,LUZ,I,FNBLNK,ISTART,ILENBAS,LINBLNK
C
      write(6,*) ' izatom on input ',izatom
C
      write(6,*) ' basnam on input '
      write(6,'(A)') basnam
      write(6,'(A80)')basnam
C
C     Construct the string.
      DO 5 I=1,80
      TITLE(I:I) = ' '
    5 CONTINUE
C
      IF(IZATOM.EQ. 1) TITLE(1:2) = 'H:'
      IF(IZATOM.EQ. 2) TITLE(1:2) = 'HE:'
      IF(IZATOM.EQ. 3) TITLE(1:2) = 'LI:'
      IF(IZATOM.EQ. 4) TITLE(1:2) = 'BE:'
      IF(IZATOM.EQ. 5) TITLE(1:2) = 'B:'
      IF(IZATOM.EQ. 6) TITLE(1:2) = 'C:'
      IF(IZATOM.EQ. 7) TITLE(1:2) = 'N:'
      IF(IZATOM.EQ. 8) TITLE(1:2) = 'O:'
      IF(IZATOM.EQ. 9) TITLE(1:2) = 'F:'
      IF(IZATOM.EQ.10) TITLE(1:2) = 'NE:'
      IF(IZATOM.EQ.11) TITLE(1:2) = 'NA:'
      IF(IZATOM.EQ.12) TITLE(1:2) = 'MG:'
      IF(IZATOM.EQ.13) TITLE(1:2) = 'AL:'
      IF(IZATOM.EQ.14) TITLE(1:2) = 'SI:'
      IF(IZATOM.EQ.15) TITLE(1:2) = 'P:'
      IF(IZATOM.EQ.16) TITLE(1:2) = 'S:'
      IF(IZATOM.EQ.17) TITLE(1:3) = 'CL:'
      IF(IZATOM.EQ.18) TITLE(1:3) = 'AR:'
C
      ISTART=linblnk(TITLE)+1
      ILENBAS=linblnk(BASNAM)
      write(6,*) ' istart, ilenbas ',istart,ilenbas

      write(6,'(A80)') title

      DO 10 I=1,ILENBAS
      TITLE(ISTART-1+I : ISTART-1+I) = BASNAM(I:I)
   10 CONTINUE
C
      WRITE(6,*) ' string to search SPLBAS for '
      write(6,'(A)') TITLE
C
C     Open the SPLBAS file.
C
      LUZ=71
      OPEN(UNIT=LUZ,FILE='SPLBAS',FORM='FORMATTED',STATUS='OLD')
C
    1 READ(LUZ,'(A)',END=5400) WRK
      IF(INDEX(WRK,':').EQ.0) GOTO 1

      write(6,*) ' lengths of characters '
      write(6,*) linblnk(wrk)
      write(6,*) linblnk(title)
      write(6,*) fnblnk(wrk)
      write(6,*) fnblnk(title)


      IF(WRK(1:LINBLNK(WRK)) .NE. TITLE(1:LINBLNK(TITLE))) GOTO 1
c     IF(WRK(1:LINBLNK(WRK)) .NE. TITLE(1:FNBLNK(TITLE)-1) ) GOTO 1
        write(6,'(A)') WRK
        write(6,'(A)') TITLE(1:ISTART+ILENBAS-1)
C
      READ(LUZ,'(A6)') BASTYP
      READ(LUZ,*) NSFULL,NPFULL,NDFULL
C
      CLOSE(UNIT=LUZ,STATUS='KEEP')
C
      write(6,*) ' bastyp '
      write(6,'(A6)') bastyp
      write(6,*) nsfull,npfull,ndfull
C
      RETURN
 5400 WRITE(6,5401)
 5401 FORMAT(T3,'@GTFLGS-F, Basis not found on SPLBAS. ')
      CLOSE(UNIT=LUZ,STATUS='KEEP')
      CALL ERREX
      END
