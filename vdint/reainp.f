      SUBROUTINE REAINP(WORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER ( NTABLE = 5 )
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
      LOGICAL BIGVC, GENCON
      COMMON /CBIREA/ IPRINT, LUMLCL, MAXPRI, BIGVC, GENCON
C
      SAVE SET
C
      DATA TABLE /'.PRINT ', '.UNIT  ', '.MAXPRI', '.BIGVEC',
     *            '.GENCON'/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize /CBIREA/
C
      IPRINT = 1
      LUMLCL = 9 
      MAXPRI=21
      BIGVC  = .FALSE.
      GENCON = .FALSE.
C
      NEWDEF = (WORD .EQ. '*READIN')
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
     *            '" not recognized in REAINP.'
               STOP
    1          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,'(I5)') LUMLCL
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,'(I5)') MAXPRI
               GO TO 100
    4          CONTINUE
                  BIGVC = .TRUE.
               GO TO 100
    5             GENCON = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in REAINP.'
               STOP 
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for READIN:',0)
         IF (IPRINT .NE. IPRDEF) THEN
            WRITE (LUPRI,'(A,I5)') ' Print level in READIN:',IPRINT
         END IF
         IF (LUMLCL .NE. LUCMD) THEN
            WRITE (LUPRI,'(/A,I5)') ' MOLECULE input from unit',LUMLCL
         END IF
         IF (MAXPRI .NE. 0) THEN
            WRITE (LUPRI,'(/A,I5)')
     *         ' Maximum number of primitives per integral block :',
     *         MAXPRI
         END IF
         IF (BIGVC) THEN
            WRITE (LUPRI,'(/A)')
     *         ' Primitives from different centers treated '//
     *         'simultaneously.'
         END IF
         IF (GENCON) THEN
            WRITE (LUPRI,'(/A)')
     *         ' Routines for general contraction used.'
         END IF
         WRITE (LUPRI,'()')
      END IF
      RETURN
      END
