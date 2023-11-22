      SUBROUTINE TWOINP(WORD)
C
C     TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (NTABLE = 17)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
      LOGICAL         SKIP,   RUNPTR, RUNSOR, RUNINT, DTEST, TKTIME,
     *                RETUR,  NODC,   NODV,   NOPV
      COMMON /CBITWO/ SKIP,   RUNPTR, RUNSOR, RUNINT, DTEST, TKTIME,
     *                IPRALL, IPRPTR, IPRSOR,
     *                IPRINT, IPRNTA, IPRNTB, IPRNTC, IPRNTD, RETUR,
     *                NODC,   NODV,   NOPV,
     *                MAXDIF, MAXVEC
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.TEST  ', '.FIRST ', '.SECOND',
     *            '.PTRSKI', '.SORSKI', '.INTSKI',
     *            '.PTRPRI', '.SORPRI', '.INTPRI',
     *            '.NODC  ', '.NODV  ', '.NOPV  ', '.RETURN',
     *            '.MAXVEC', '.TIME  '/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize /CBITWO/
C
      SKIP   = .FALSE.
      RUNPTR = .TRUE.
      RUNSOR = .TRUE.
      RUNINT = .TRUE.
      IPRALL = IPRDEF
      IF (MOLHES) THEN
         MAXDIF = 2
      ELSE IF (MOLGRD) THEN
         MAXDIF = 1
      ELSE
         SKIP = .TRUE.
      END IF
      IPRINT = IPRDEF
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      IPRPTR = IPRDEF
      IPRSOR = IPRDEF
      NODC   = .FALSE.
      NODV   = .FALSE.
      NOPV   = .FALSE.
      DTEST  = .FALSE.
      RETUR  = .FALSE.
      MAXVEC = 8
      TKTIME = .FALSE.
C
      NEWDEF = WORD .EQ. '*TWOEXP'
      MAXOLD = MAXDIF
      ICHANG = 0
      IF (NEWDEF) THEN
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),
     *                      I
                  END IF
  200          CONTINUE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in TWOINP.'
               STOP
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRALL
                  IF (IPRALL .EQ. IPRDEF) ICHANG = ICHANG - 1
                  IPRINT = IPRALL
                  IPRPTR = IPRALL
                  IPRSOR = IPRALL
               GO TO 100
    3             DTEST = .TRUE.
               GO TO 100
    4             MAXDIF = 1
                  IF (MAXDIF .EQ. MAXOLD) ICHANG = ICHANG - 1
               GO TO 100
    5             MAXDIF = 2
                  IF (MAXDIF .EQ. MAXOLD) ICHANG = ICHANG - 1
               GO TO 100
    6             RUNPTR = .FALSE.
               GO TO 100
    7             RUNSOR = .FALSE.
               GO TO 100
    8             RUNINT = .FALSE.
               GO TO 100
    9             READ (LUCMD, '(I5)') IPRPTR
                  IF (IPRPTR .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   10             READ (LUCMD, '(I5)') IPRSOR
                  IF (IPRSOR .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   11             READ (LUCMD, '(5I5)') IPRINT, IPRNTA, IPRNTB,
     *                                          IPRNTC, IPRNTD
                  IPRSUM = IPRNTA + IPRNTB + IPRNTC + IPRNTD
                  IF (IPRINT .EQ. IPRDEF .AND. IPRSUM .EQ. 0) THEN
                     ICHANG = ICHANG - 1
                  END IF
               GO TO 100
   12             NODC = .TRUE.
               GO TO 100
   13             NODV = .TRUE.
               GO TO 100
   14             NOPV = .TRUE.
               GO TO 100
   15             RETUR = .TRUE.
               GO TO 100
   16             READ (LUCMD, '(I5)') MAXVEC
                  WRITE (LUPRI,'(/A/)')
     *             ' .MAXVEC not implemented in TWOINP - '//
     *             'use .MAXPRI under *READIN instead.'
                  STOP
C              GO TO 100
   17             TKTIME = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in TWOINP.'
               STOP 
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .EQ. 0) RETURN
      IF (NEWDEF) THEN
         CALL HEADER('Changes of defaults for TWOEXP:',0)
      END IF
      IF (SKIP) THEN
         WRITE (LUPRI,'(A)') ' TWOEXP skipped in this run.'
      ELSE IF (NEWDEF) THEN
         IF (IPRALL .NE. IPRDEF) THEN
            WRITE (LUPRI,'(A,I5)') ' Print level in TWOEXP:',IPRALL
         END IF
         IF (RUNPTR) THEN
            IF (IPRPTR .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in PTRAN :',IPRPTR
            END IF
         ELSE
            WRITE (LUPRI,'(A)') ' PTRAN skipped in this run.'
         END IF
         IF (RUNSOR) THEN
            IF (IPRSOR .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in PSORT :',IPRSOR
            END IF
         ELSE
            WRITE (LUPRI,'(A)') ' PSORT skipped in this run.'
         END IF
         IF (RUNINT) THEN
            IF (MAXDIF .NE. MAXOLD) THEN
               WRITE (LUPRI,'(A,I1)') ' Maximum differentiation: ',
     *                                MAXDIF
            END IF
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in TWOINT:',
     *                                IPRINT
            END IF
            IF (IPRNTA + IPRNTB + IPRNTC + IPRNTD .GT. 0) THEN
               WRITE(LUPRI,'(2A,4I3)')' Extra output for the following',
     *                    ' shells:', IPRNTA, IPRNTB, IPRNTC, IPRNTD
               IF (RETUR) WRITE (LUPRI,'(A)')
     *             ' Program will exit TWOINT after these shells.'
            END IF
            IF (MAXVEC .NE. 8) THEN
               WRITE(LUPRI,'(A,I5)') ' MAXVEC has been set to ',MAXVEC
            END IF
            IF (DTEST) THEN
               WRITE (LUPRI,'(/,2A)')
     *           ' Test run for two-electron integrals',
     *           '- special input required by DTEST0 call from TWOINT.'
            END IF
            IF (NODC) WRITE (LUPRI,'(/,2A)') ' Inactive one-electron',
     *         ' density matrix neglected in TWOEXP.'
            IF (NODV) WRITE (LUPRI,'(/,2A)') ' Active one-electron',
     *         ' density matrix neglected in TWOEXP.'
            IF (NOPV) WRITE (LUPRI,'(/,2A)') ' Active two-electron',
     *         ' density matrix neglected in TWOEXP.'
            IF (TKTIME) WRITE (LUPRI,'(/,2A)') ' Detailed timing for',
     *         ' integral calculation will be provided.'
         ELSE
            WRITE (LUPRI,'(A)') ' TWOINT skipped in this run.'
         END IF
      END IF
      RETURN
      END
