      SUBROUTINE DIPINP(WORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (NTABLE = 5)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
      LOGICAL SKIP, NODC, NODV, TEST
      COMMON /CBIDIP/ IPRINT, SKIP, NODC, NODV, TEST
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.NODC  ', '.NODV  ',
     *            '.TEST  '/
      DATA SET/.FALSE./
C
      SKIP = .NOT.(DIPDER .OR. POLAR)
      IF (SET) RETURN
C     IF (SET) THEN
C        IF (.NOT. SKIP) SKIP = .NOT.(DIPDER .OR. POLAR)
C        RETURN
C     END IF
      SET = .TRUE.
C
C     Initialize /CBIDIP/
C
      IPRINT = IPRDEF
      NODC   = .FALSE.
      NODV   = .FALSE.
      TEST   = .FALSE.
C
      NEWDEF = (WORD .EQ. '*DIPCTL')
      ICHANG = 0
      IF (NEWDEF) THEN
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5), I
                  END IF
  200          CONTINUE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in DIPINP.'
               STOP 
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  NODC = .TRUE.
               GO TO 100
    4          CONTINUE
                  NODV = .TRUE.
               GO TO 100
    5             TEST = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in DIPINP.'
               STOP 
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for DIPCTL:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' DIPCTL skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in DIPCTL:',IPRINT
            END IF
            IF (NODC) WRITE (LUPRI,'(/,2A)') ' Inactive one-electron',
     *      ' density matrix neglected in DIPCTL.'
            IF (NODV) WRITE (LUPRI,'(/,2A)') ' Active one-electron',
     *      ' density matrix neglected in DIPCTL.'
            IF (TEST) WRITE (LUPRI,'(/,2A)') ' Test for dipole moments',
     *      ' and dipole reorthonormalization.'
         END IF
      END IF
      RETURN
      END
