      SUBROUTINE ORTINP(WORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (NTABLE = 2)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
      LOGICAL SKIP
      COMMON /CBIORT/ IPRINT, SKIP
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT '/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize /CBIORT/
C
      IPRINT = IPRDEF
      SKIP   = .NOT.MOLHES
C
      NEWDEF = (WORD .EQ. '*REORT ')
      ICHANG = 0
      IF (NEWDEF) THEN
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2), I
                  END IF
  200          CONTINUE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ORTINP.'
               STOP
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ORTINP.'
               STOP 
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for REORT:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' REORT skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in REORT :',IPRINT
            END IF
         END IF
      END IF
      RETURN
      END
