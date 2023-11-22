      SUBROUTINE RESINP(WORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
      PARAMETER (NTABLE = 13)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      LOGICAL SKIP, NOTRIA, NEWRD, DOXKAP, RSTNR, NOAVER, D1DIAG, PRIL2
      COMMON /CBIRES/ THRRNR, MAXRNR,
     *                IPRINT, MAXSIM, IRDC, NRDC, NRDT, NRDCO(MXCOOR),
     *                SKIP, NOTRIA, NEWRD, DOXKAP, RSTNR, NOAVER,
     *                D1DIAG, PRIL2
C
C     Used from common blocks:
C       NRINF  : THRNR, MAXNR
C       NUCLEI : NDCORD
C
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
C     NRINF : space nedded in NRLIN = LNRWA + nsim*LNRWB
C             NSIDE=0, symmetric; =1 from left side; =2 from right side
C             NREFS = # of reference vectors to othogonalize against
      COMMON /NRINF / THRNR, MAXNR, IPRNR, NRTYPA,NRTYPB,NRTOT, NSIDE,
     *                NREFS, LUNRA, LUNRB, LUNRC,
     *                LNRWA, LNRWB, KNRLIN(20),   NNRLIN,LNRLIN
      CHARACTER NAMEX*6
      LOGICAL DCORD, DCORGD, NOORBT, DOPERT
      COMMON /NUCLEIi/ NOORBT(MXCENT),
     &                NUCIND, NUCDEP, NUCPRE(MXCENT), NUCNUM(MXCENT,8),
     &                NUCDEG(MXCENT), ISTBNU(MXCENT), NDCORD,
     &                NDCOOR(MXCOOR), NTRACO, NROTCO, ITRACO(3),
     &                IROTCO(3),
     &                NATOMS, NFLOAT,
     &                IPTGDV(3*MXCENT),
     &                NGDVEC(8), IGDVEC(8)
      COMMON /NUCLEI/ CHARGE(MXCENT), CORD(MXCENT,3),
     &                DCORD(MXCENT,3),DCORGD(MXCENT,3),
     &                DOPERT(0:3*MXCENT)
      COMMON /NUCLEC/ NAMEX(MXCOOR)
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ','.MAXITE','.THRESH',
     *            '.MAXSIM', '.NEWRD ','.NOTRIA','.RDVECS',
     *            '.NONEXT', '.NRREST','.NOAVER','.D1DIAG',
     *            '.PRIL2 '/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize /CBIRES/
C
      IPRINT = IPRDEF
      MAXSIM = 15
      NRDT   = 0
      SKIP   = .FALSE.
      NOTRIA = .FALSE.
      NEWRD  = .FALSE.
      DOXKAP = .TRUE.
      RSTNR  = .FALSE.
      NOAVER = .FALSE.
      D1DIAG = .FALSE.
      PRIL2  = .FALSE.
C
C     Initialize /NRINF/
C
      THRNR  = 1.D-03
      MAXNR  = 20
C
      NEWDEF = (WORD .EQ. '*RESPON')
      THRNRD = THRNR
      MAXNRD = MAXNR
      MAXSMD = MAXSIM
      IF (.NOT. (MOLHES .OR. DIPDER .OR. POLAR)) SKIP = .TRUE.
      ICHANG = 0
      IF (NEWDEF) THEN
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13), I
                  END IF
  200          CONTINUE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in RESINP.'
               STOP
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,'(I5)') MAXNR
                  IF (MAXNR .EQ. MAXNRD) ICHANG = ICHANG - 1
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,'(E12.6)') THRNR
                  IF (THRNR .EQ. THRNRD) ICHANG = ICHANG - 1
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,'(I5)') MAXSIM
                  IF (MAXSIM .EQ. MAXSMD)ICHANG = ICHANG - 1
               GO TO 100
    6          CONTINUE
                  NEWRD  = .TRUE.
               GO TO 100
    7             NOTRIA = .TRUE.
               GO TO 100
    8          CONTINUE
                  READ (LUCMD,*) NRDT
                  READ (LUCMD,*) (NRDCO(I), I = 1,NRDT)
               GO TO 100
    9             DOXKAP = .FALSE.
               GO TO 100
   10             RSTNR = .TRUE.
               GO TO 100
   11             NOAVER = .TRUE.
               GO TO 100
   12             D1DIAG = .TRUE.
               GO TO 100
   13             PRIL2  = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in RESINP.'
               STOP
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for RESPON:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' RESPON skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     *         ' Print level in RESPON        :',IPRINT
            IF (THRNR .NE. THRNRD) WRITE (LUPRI,'(A,1P,E9.2)')
     *            ' Threshold in RESPON          :',THRNR
            IF (MAXNRD .NE. MAXNR) THEN
               WRITE(LUPRI,'(A,I5)')' Maximum iterations in RESPON :',
     *                              MAXNR
            END IF
            IF (MAXSIM .NE. MAXSMD) THEN
               WRITE(LUPRI,'(A,I5)')' Maximum vectors in NRCTL     :',
     *                              MAXSIM
            END IF
            IF (NEWRD) THEN
               WRITE (LUPRI,'(/A)') ' New RD file will be used.'
            ELSE IF (NOTRIA) THEN
               WRITE (LUPRI,'(/2A)') ' Old solution vectors not',
     *            ' used as trial vectors in RESPON.'
            END IF
            IF (NRDT .GT. 0) THEN
               WRITE (LUPRI,'(/A,(T45,10I3))')
     *            ' Only solve specified gradient vectors:',
     *            (NRDCO(I), I = 1,NRDT)
            END IF
            IF (.NOT.DOXKAP) WRITE (LUPRI,'(/A)')
     *         ' Straightforward Conjugate Gradient in this run!'
            IF (RSTNR) WRITE (LUPRI,'(/A)')
     *         ' Restart of RESPON from saved trial vectors.'
            IF (D1DIAG) WRITE (LUPRI,'(/2A)')
     *         ' Diagonal Hessian elements not used when',
     *         ' generating trial vectors.'
            IF (NOAVER) WRITE (LUPRI,'(/2A)')
     *         ' Diagonal orbital Hessian for trial vectors',
     *         ' approximated by Fock contributions.'
            IF (PRIL2) WRITE (LUPRI,'(/2A)')
     *         ' The electronic Hessian will be calculated explicitly',
     *         ' and tested for symmetry.'
         END IF
      END IF
      IF (NRDT .LT. 0 .OR. NRDT .GT. NDCORD+6) THEN
         WRITE (LUPRI,'(A,I10)') ' *** RESINP ERROR: Illegal NRDT:',NRDT
         STOP
      END IF
      IPRNR = IPRINT
      RETURN
      END
