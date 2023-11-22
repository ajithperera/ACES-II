      SUBROUTINE CRSINP(WORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
      LOGICAL SET, NEWDEF
      PARAMETER (NTABLE = 41)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      DIMENSION NATOM(10)
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
      LOGICAL SKIP, RUNINT,
     *        RUNFCK, TWOH1, NOH1, NOH2, FCKTST, RUNSDR, RUNTRA,
     *        SDRTST, UNDIF, TRATST, RUNGDY, NODC, NODV, NOPV, TESTFS,
     *        DIAGTD, NOFD , NOFS, NOSSF, NOORTH, J1TST, RETUR, TKTIME,
     *        AOMAT, NOKAP1, OITST, RUNOIT, RHSTES, RUNGJ1
      COMMON /CBICRS/ SKIP, RUNINT, RUNFCK, IPRALL, IPRFCK, IPRINT,
     *                IPRNTA, IPRNTB, IPRNTC, IPRNTD, IPRMAX, MAXSIM,
     *                TWOH1, NOH1, NOH2, FCKTST, RUNSDR, RUNTRA,
     *                IPRSDR, IPRTRA, SDRTST, UNDIF, TRATST, RUNGDY,
     *                IPRGDY, NODC, NODV, NOPV, TESTFS, DIAGTD, NOFD,
     *                NOFS, NOSSF, NOORTH, J1TST, RETUR, MAXVEC, TKTIME,
     *                AOMAT, NOKAP1, OITST, RUNOIT, RHSTES, RUNGJ1
      LOGICAL CCTFLG
      COMMON /CCTGRD/ CCTDOT(MXCOOR), CCTEXP(MXCOOR), CCTFLG(MXCOOR)
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ', '.MAXSIM',
     *            '.INTSKI', '.INTPRI', '.RETURN',
     *            '.FCKSKI', '.FCKPRI', '.TWOH1 ',
     *            '.NOH1  ', '.NOH2  ', '.FCKTES',
     *            '.SDRSKI', '.TRASKI', '.SDRPRI',
     *            '.TRAPRI', '.SDRTES', '.UNDIF ',
     *            '.TRATES', '.GDYSKI', '.GDYPRI',
     *            '.NODC  ', '.NODV  ', '.NOPV  ',
     *            '.FSTTES', '.DIAGTD', '.NOFD  ',
     *            '.NOFS  ', '.NOSSF ', '.SIRPR4',
     *            '.SIRPR6', '.NOORTH', '.J1TEST',
     *            '.MAXVEC', '.TIME  ', '.AOMAT ',
     *            '.NOKAP1', '.OITEST', '.RHSSKI',
     *            '.RHSTES', '.GJ1SKI'/
      DATA SET/.FALSE./, IPRI4/0/, IPRI6/0/
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize /CBICRS/
C
      SKIP   = .NOT. (MOLHES .OR. DIPDER)
      RUNINT = .TRUE.
      RUNFCK = .TRUE.
      IPRALL = IPRDEF
      IPRFCK = IPRDEF
      IPRINT = IPRDEF
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      MAXSIM = 3
      TWOH1  = .FALSE.
      NOH1   = .FALSE.
      NOH2   = .FALSE.
      FCKTST = .FALSE.
      RUNSDR = .TRUE.
      RUNTRA = .TRUE.
      RUNOIT = .TRUE.
      RUNGJ1 = .TRUE.
      IPRSDR = IPRDEF
      IPRTRA = IPRDEF
      SDRTST = .FALSE.
      UNDIF  = .FALSE.
      TRATST = .FALSE.
      RUNGDY = MOLHES
      IPRGDY = IPRDEF
      NODC   = .FALSE.
      NODV   = .FALSE.
      NOPV   = .FALSE.
      TESTFS = .FALSE.
      RHSTES = .FALSE.
      DIAGTD = .FALSE.
      NOFD   = .FALSE.
      NOFS   = .FALSE.
      NOSSF  = .FALSE.
      NOORTH = .FALSE.
      OITST  = .FALSE.
      RETUR  = .FALSE.
      MAXVEC = 8
      TKTIME = .FALSE.
      AOMAT  = .FALSE.
      NOKAP1 = .FALSE.
C
      NEWDEF = WORD .EQ. '*CCRHS '
      ICHANG = 0
      IF (NEWDEF) THEN
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     ICHANG = ICHANG + 1
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,
     *                      15,16,17,18,19,20,21,22,23,24,25,
     *                      26,27,28,29,30,31,32,33,34,35,36,
     *                      37,38,39,40,41), I
                  END IF
  200          CONTINUE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in CRSINP.'
               STOP 
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRALL
                  IF (IPRALL .EQ. IPRDEF) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     IF (IPRINT .EQ. IPRDEF) IPRINT = IPRALL
                     IF (IPRFCK .EQ. IPRDEF) IPRFCK = IPRALL
                     IF (IPRSDR .EQ. IPRDEF) IPRSDR = IPRALL
                     IF (IPRTRA .EQ. IPRDEF) IPRTRA = IPRALL
                     IF (IPRGDY .EQ. IPRDEF) IPRGDY = IPRALL
                  END IF
               GO TO 100
    3             READ (LUCMD, '(I5)') MAXSIM
                  IF (MAXSIM .EQ. 3) ICHANG = ICHANG - 1
               GO TO 100
    4             RUNINT = .FALSE.
               GO TO 100
    5             READ (LUCMD, '(5I5)') IPRINT, IPRNTA, IPRNTB,
     *                                      IPRNTC, IPRNTD
                  IPRSUM = IPRNTA + IPRNTB + IPRNTC + IPRNTD
                  IF (IPRINT .EQ. IPRDEF .AND. IPRSUM .EQ. 0) THEN
                     ICHANG = ICHANG - 1
                  END IF
               GO TO 100
    6             CONTINUE
                  RETUR = .TRUE.
               GO TO 100
    7             RUNFCK = .FALSE.
               GO TO 100
    8             READ (LUCMD, '(I5)') IPRFCK
                  IF (IPRFCK .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    9             TWOH1 = .TRUE.
               GO TO 100
   10             NOH1 = .TRUE.
               GO TO 100
   11             NOH2 = .TRUE.
               GO TO 100
   12             FCKTST = .TRUE.
               GO TO 100
   13             RUNSDR = .FALSE.
               GO TO 100
   14             RUNTRA = .FALSE.
               GO TO 100
   15             READ (LUCMD, '(I5)') IPRSDR
                  IF (IPRSDR .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   16             READ (LUCMD, '(I5)') IPRTRA
                  IF (IPRTRA .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   17             SDRTST = .TRUE.
               GO TO 100
   18             UNDIF = .TRUE.
               GO TO 100
   19             TRATST = .TRUE.
               GO TO 100
   20             RUNGDY = .FALSE.
               GO TO 100
   21             READ (LUCMD, '(I5)') IPRGDY
                  IF (IPRGDY .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   22             NODC = .TRUE.
               GO TO 100
   23             NODV = .TRUE.
               GO TO 100
   24             NOPV = .TRUE.
               GO TO 100
   25             TESTFS = .TRUE.
               GO TO 100
   26             DIAGTD = .TRUE.
               GO TO 100
   27             NOFD = .TRUE.
               GO TO 100
   28             NOFS = .TRUE.
               GO TO 100
   29             NOSSF = .TRUE.
               GO TO 100
   30             READ (LUCMD,*) IPRI4
                  IF (IPRI4 .EQ. 0) ICHANG = ICHANG - 1
               GO TO 100
   31             READ (LUCMD,*) IPRI6
                  IF (IPRI6 .EQ. 0) ICHANG = ICHANG - 1
               GO TO 100
   32             NOORTH = .TRUE.
               GO TO 100
   33             J1TST  = .TRUE.
               GO TO 100
   34             READ (LUCMD,'(I5)') MAXVEC
                  WRITE (LUPRI,'(/A/)')
     *             ' .MAXVEC not implemented in CRSINP - '//
     *             'use .MAXPRI under *READIN instead.'
                  STOP 
C              GO TO 100
   35             TKTIME  = .TRUE.
               GO TO 100
   36             AOMAT   = .TRUE.
               GO TO 100
   37             NOKAP1  = .TRUE.
               GO TO 100
   38             OITST   = .TRUE.
               GO TO 100
   39             RUNOIT  = .FALSE.
               GO TO 100
   40             RHSTES  = .TRUE.
               GO TO 100
   41             RUNGJ1  = .FALSE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in CRSINP.'
               STOP 
            END IF
      END IF
  300 CONTINUE
C
      IPRMAX = MAX(IPRALL,IPRINT,IPRFCK,IPRSDR,IPRTRA,IPRGDY)
C
      IF (NEWDEF .AND. (ICHANG .GT. 0) .AND. SKIP) THEN
         CALL HEADER('Changes of defaults for CRSIDE:',0)
         WRITE (LUPRI,'(A)') ' CRSIDE skipped in this run.'
      ELSE IF (NEWDEF .AND. ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for CRSIDE:',0)
         IF (IPRALL .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     *         ' Print level in CRSIDE:',IPRALL
         IF (TKTIME) WRITE (LUPRI,'(A)')
     *         ' Detailed timing for integral calculation '//
     *          'will be provided.'
         IF (NODC) WRITE (LUPRI,'(A)')
     *         ' DC will be neglected in this run.'
         IF (NODV) WRITE (LUPRI,'(A)')
     *         ' DV will be neglected in this run.'
         IF (NOPV) WRITE (LUPRI,'(A)')
     *         ' PV will be neglected in this run.'
         IF (UNDIF) WRITE (LUPRI,'(A)')
     *         ' Undifferentiated integrals replaces J1.'
         IF (NOH1) WRITE (LUPRI,'(A)')
     *         ' One-electron Hamiltonian neglected.'
         IF (NOH2) WRITE (LUPRI,'(A)')
     *         ' Two-electron Hamiltonian neglected.'
         IF (NOORTH) WRITE (LUPRI,'(A)')
     *         ' No orthonormalization in differentiated gradients.'
         IF (J1TST)  WRITE (LUPRI,'(A)')
     *         ' Test performed for J1 integrals.'
         IF (OITST)  WRITE (LUPRI,'(A)')
     *         ' Test performed for one-index transformations.'
         IF (.NOT.RUNOIT) WRITE (LUPRI,'(A)')
     *         ' No right-hand side constructed (i.e., skipping ONEIT).'
         IF (.NOT.RUNGJ1) WRITE (LUPRI,'(A)') ' CCGDJ1 skipped.'
         IF (RHSTES) WRITE (LUPRI,'(A)')
     *         ' Test performed for right-hand side.'
         IF (NOKAP1)  WRITE (LUPRI,'(A)')
     *         ' Kappa(1) set to zero.'
         IF (MAXSIM .NE. 3) WRITE (LUPRI,'(A,I3)')
     *         ' Maximum number of simultaneous directions:',MAXSIM
C
C        TWOINT
C
         IF (RUNINT) THEN
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in TWOINT:',IPRINT
            END IF
            IF (IPRNTA + IPRNTB + IPRNTC + IPRNTD .GT. 0) THEN
               WRITE (LUPRI,'(2A,4I3)') ' Extra output for the',
     *            ' following shells:', IPRNTA, IPRNTB, IPRNTC, IPRNTD
               IF (RETUR) WRITE (LUPRI,'(A)')
     *            'Program will exit TWOINT after these shells.'
            END IF
            IF (MAXVEC .NE. 8) THEN
               WRITE (LUPRI,'(A,I5)') ' MAXVEC has been set to ',MAXVEC
            END IF
         ELSE
            WRITE (LUPRI,'(A)') ' TWOINT skipped in this run.'
         END IF
C
C        GETFD
C
         IF (RUNFCK) THEN
            IF (IPRFCK .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in GETFD :',IPRFCK
            END IF
            IF (FCKTST) WRITE (LUPRI, '(A)')
     *         ' Testing of Fock matrices.'
            IF (AOMAT) WRITE (LUPRI, '(A)')
     *         ' AO printing of SO Fock matrices.'
            IF (TWOH1) THEN
               WRITE (LUPRI, '(A)')
     *           ' One-electron Hamiltonian multiplied by two in GETFD.'
               WRITE (LUPRI, '(2A)') ' WARNING: TWOH1 not implemented',
     *            ' in this version.'
               STOP 
            END IF
         ELSE
            WRITE (LUPRI,'(A)') ' GETFD skipped in this run.'
         END IF
C
C        GETSD
C
         IF (RUNSDR) THEN
            IF (IPRSDR .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in GETSD :',IPRSDR
            END IF
            IF (DIAGTD) WRITE (LUPRI, '(A)')
     *          ' Diagonal overlap matrices (TD(I,I) = 0.25)'
            IF (SDRTST) WRITE (LUPRI, '(A)')
     *         ' AO printing of SO overlap matrices.'
            IF (AOMAT) WRITE (LUPRI, '(A)')
     *       ' SO overlap matrices transformed and printed in AO basis.'
         ELSE
            WRITE (LUPRI,'(A)') ' GETSD skipped in this run.'
         END IF
C
C        DERTRA
C
         IF (RUNTRA) THEN
            IF (IPRTRA .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in DERTRA:',IPRTRA
            END IF
            IF (TRATST) WRITE (LUPRI,'(A)')' Testing on transformation.'
         ELSE
            WRITE (LUPRI,'(A)') ' DERTRA skipped in this run.'
         END IF
C
C        GETGDY
C
         IF (RUNGDY) THEN
            IF (IPRGDY .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in GETGDY:',IPRGDY
            END IF
            IF (TESTFS) WRITE (LUPRI,'(A)')
     *         ' Testing on one-index transformation of Fock matrices'
            IF (NOFD) WRITE (LUPRI,'(A)')
     *         ' FD not included in calculation of Y matrices'
            IF (NOFS) WRITE (LUPRI,'(A)')
     *         ' FS not included in calculation of Y matrices'
            IF (NOSSF) WRITE (LUPRI,'(A)')
     *         ' SSF not included in calculation of Y matrices'
         ELSE
            WRITE (LUPRI,'(A)') ' GETGDY skipped in this run.'
         END IF
      END IF
C
C     Initialize /CCTGRD/
C
      DO 400 ICOOR = 1, MXCOOR
         CCTFLG(ICOOR) = .FALSE.
  400 CONTINUE
      RETURN
      END
