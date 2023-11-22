      SUBROUTINE ONEINP(WORD)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (NTABLE = 5)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     &        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     &        H2MO
      COMMON /ABAINF/ IPRDEF,
     &                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     &                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     &                H2MO
      LOGICAL         SKIP, DTEST, DIFINT, NODC, NODV, DIFDIP
      COMMON /CBIONE/ IPRINT, SKIP, MAXDIF, DTEST, IDATOM,
     &                IDCOOR, DIFINT, NODC, NODV, DIFDIP
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ','.TEST  ',
     &            '.NODC  ', '.NODV  '/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C INITIALIZE /CBIONE/
C
      IPRINT = IPRDEF
      SKIP   = .FALSE.
      IF (MOLHES .OR. DIPDER) THEN
         MAXDIF = 2
      ELSE IF (MOLGRD .OR. POLAR) THEN
         MAXDIF = 1
      ELSE
         SKIP = .TRUE.
      END IF
      DTEST  = .FALSE.
      NODC   = .FALSE.
      NODV   = .FALSE.
      DIFDIP = DIPDER .OR. POLAR
C
      NEWDEF = WORD .EQ. '*ONEINT'
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
     *            '" not recognized in ONEINP.'
               STOP
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3             DTEST = .TRUE.
                  READ (LUCMD, '(2I5,L5)') IDATOM, IDCOOR, DIFINT
               GO TO 100
    4             NODC = .TRUE.
               GO TO 100
    5             NODV = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ONEINP.'
               STOP
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for ONEINT:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' ONEINT skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in ONEINT:',IPRINT
            END IF
            IF (NODC) WRITE (LUPRI,'(/,2A)') ' Inactive one-electron',
     *      ' density matrix neglected in ONEINT.'
            IF (NODV) WRITE (LUPRI,'(/,2A)') ' Active one-electron',
     *      ' density matrix neglected in ONEINT.'
            IF (DTEST) THEN
               IF (.NOT. DIFINT) THEN
                  WRITE (LUPRI,'(/,2A)') ' Numerical testing',
     *                        ' of one-electron expectation values'
               ELSE
                  WRITE (LUPRI,'(/,2A)') ' Numerical testing of',
     *              ' one-electron expectation values and integrals'
               END IF
               WRITE(LUPRI,'(/,A,2I3)') ' Atom and coordinate: ',
     *                                    IDATOM, IDCOOR
            END IF
         END IF
      END IF
      RETURN
      END
