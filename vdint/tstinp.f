
      SUBROUTINE TSTINP(WORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SET, SKIP, NEWDEF, INTTST, RETUR
      PARAMETER (NTABLE = 9)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
      COMMON /CBITST/ SKIP, IPRINT, IPRNTA, IPRNTB, IPRNTC, IPRNTD,
     *                IATOM, ICOOR, DELTA, TSTTHR, NORDER, INTTST,
     *                RETUR, MAXVEC
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.PERTUR',
     *            '.DELTA ', '.THRESH', '.SECOND','.INTEGR',
     *            '.RETURN','.MAXVEC'/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize CBITST
C
      SKIP   = .TRUE.
      IPRINT = IPRDEF
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      DELTA  = 1.0D-05
      TSTTHR = 1.0D-10
      NORDER = 1
      INTTST = .FALSE.
      RETUR  = .FALSE.
      MAXVEC = 8
C
      NEWDEF = WORD .EQ. '*TWOTST'
      IF (NEWDEF) THEN
         SKIP = .FALSE.
         CALL HEADER('CHANGES OF DEFAULTS FOR TWOTST:',0)
      END IF
      IF (NEWDEF) THEN
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9), I
                  END IF
  200          CONTINUE
               WRITE (LUPRI,'(/,3A,/)') ' KEYWORD "',WORD,
     *            '" NOT RECOGNIZED IN TSTINP.'
               STOP
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(5I5)') IPRINT, IPRNTA, IPRNTB,
     *                                          IPRNTC, IPRNTD
               GO TO 100
    3          CONTINUE
                  READ (LUCMD, '(2I5)') IATOM, ICOOR
               GO TO 100
    4             READ (LUCMD, '(E12.6)') DELTA
               GO TO 100
    5             READ (LUCMD, '(E12.6)') TSTTHR
               GO TO 100
    6             NORDER = 2
                  DELTA = 1.0D-04
               GO TO 100
    7             INTTST = .TRUE.
               GO TO 100
    8             RETUR = .TRUE.
               GO TO 100
    9             READ (LUCMD, '(I5)') MAXVEC
                  WRITE (LUPRI,'(/A/)')
     *             ' .MAXVEC not implemented in TSTINP - '//
     *             'use .MAXPRI under *READIN instead.'
                  STOP
C              GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' PROMPT "',WORD,
     *            '" NOT RECOGNIZED IN TSTINP.'
               STOP
            END IF
  300    CONTINUE
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' TWOTST SKIPPED IN THIS RUN.'
            RETURN
         END IF
         IF (INTTST) THEN
            WRITE (LUPRI,'(A)') ' NUMERICAL TESTING OF INTEGRALS'
         ELSE
            WRITE(LUPRI,'(A)')' NUMERICAL TESTING OF EXPECTATION VALUES'
         END IF
         WRITE (LUPRI,'(A,2(A,I2))') ' NUMERICAL DIFFERENTIATION FOR',
     *      ' ATOM',IATOM,' AND COORDINATE',ICOOR
         IF (INTTST .AND. NORDER .EQ. 2) THEN
            WRITE (LUPRI,'(A)')
     *         ' SECOND ORDER ACCURACY IN NUMERICAL DIFFERENTIATION.'
         END IF
         WRITE (LUPRI,'(A,I5)') ' PRINT LEVEL IN TWOTST:', IPRINT
         IF (IPRNTA + IPRNTB + IPRNTC + IPRNTD .GT. 0) THEN
               WRITE(LUPRI,'(2A,4I3)')' EXTRA OUTPUT FOR THE FOLLOWING',
     *           ' SHELLS:', IPRNTA, IPRNTB, IPRNTC, IPRNTD
            IF (RETUR) WRITE (LUPRI,'(A)')
     *         ' Program will exit from TWOINT after these shells.'
         END IF
         WRITE (LUPRI,'(A,E12.6)') ' DELTA IS SET TO     ',DELTA
         WRITE (LUPRI,'(A,E12.6)') ' THRESHOLD IS SET TO ',TSTTHR
         IF (MAXVEC .NE. 8) THEN
            WRITE (LUPRI,'(A,I5)') ' MAXVEC has been set to ',MAXVEC
         END IF
      END IF
      RETURN
      END
