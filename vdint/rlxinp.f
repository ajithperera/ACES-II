      SUBROUTINE RLXINP(WORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (NTABLE = 4)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
      LOGICAL SKIP, SYMTST, DOSELL
      COMMON /CBIRLX/ IPRINT, SKIP, SYMTST, DOSELL
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ','.SYMTES','.NOSELL'/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize /CBIRLX/
C
      IPRINT = IPRDEF
      SKIP   = .NOT. MOLHES
      SYMTST = .FALSE.
      DOSELL = .TRUE.
C
      NEWDEF = (WORD .EQ. '*RELAX')
      ICHANG = 0
      IF (NEWDEF) THEN
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4), I
                  END IF
  200          CONTINUE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in RLXINP.'
               STOP 
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  SYMTST = .TRUE.
               GO TO 100
    4             DOSELL = .FALSE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in RLXINP.'
               STOP
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for RELAX:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' RELAX skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in RELAX :',IPRINT
            END IF
            IF (SYMTST) WRITE (LUPRI,'(2A)')
     *                    ' Calculation of all elements',
     *                    ' (i,j) and (j,i) in relaxation Hessian '
            IF (.NOT.DOSELL) WRITE (LUPRI,'(A,/,2A,/)')
     *           ' Sellers'' method for quadratic errors not used.',
     *           ' Errors in relaxation may be linear in threshold',
     *           ' for response equations.'
         END IF
      END IF
      RETURN
      END
