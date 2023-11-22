      SUBROUTINE CZRINP(WORD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
      PARAMETER (NTABLE = 13)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7
      LOGICAL SKIP, EXPJAC, LRTEST, NEWRD, NOTRIA, NOSOL, RDTEST
      COMMON /CBICZRi/ SKIP, IPRINT, MAXSIM, EXPJAC, LRTEST,
     *                 NEWRD, NOTRIA, MAXIT, NOSOL, RDTEST
      COMMON /CBICZR/ DTEST, THRESH
C
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
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
     *            '.XXXXXX', '.NEWRD ','.NOTRIA','.JACOBI',
     *            '.LRSOLV', '.NOSOL ','.RDTEST','.XXXXXX',
     *            '.XXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) RETURN
      SET = .TRUE.
C
C     Initialize /CBICZR/
C
      IPRINT = IPRDEF
      MAXSIM = 1
      MAXIT  = 20
      DTEST  = 1.D-04
      SKIP   = .NOT.MOLHES
      RDTEST = .FALSE.
      NOTRIA = .FALSE.
      NEWRD  = .FALSE.
      THRESH = 1.D-06
      EXPJAC = .FALSE.
      LRTEST = .FALSE.
      NOSOL  = .FALSE.
C
      NEWDEF = (WORD .EQ. '*CCZRES')
      THRNRD = THRESH
      MAXNRD = MAXIT
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
     *            '" not recognized in CZRINP.'
               STOP 
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,'(I5)') MAXIT
                  IF (MAXIT .EQ. MAXNRD) ICHANG = ICHANG - 1
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,'(E12.6)') THRESH
                  IF (THRESH .EQ. THRNRD) ICHANG = ICHANG - 1
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
                  EXPJAC = .TRUE.
               GO TO 100
    9          CONTINUE
                  LRTEST = .TRUE.
               GO TO 100
   10          CONTINUE
                  NOSOL  = .TRUE.
               GO TO 100
   11          CONTINUE
                  RDTEST = .TRUE.
               GO TO 100
   12          CONTINUE
               GO TO 100
   13          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in CZRINP.'
               STOP 
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for CCZRES:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' CCZRES skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     *         ' Print level in CCZRES        :',IPRINT
            IF (THRESH .NE. THRNRD) WRITE (LUPRI,'(A,1P,E9.2)')
     *            ' Threshold in CCZRES          :',THRESH
            IF (MAXNRD .NE. MAXIT) THEN
               WRITE(LUPRI,'(A,I5)')' Maximum iterations in CCZRES :',
     *                              MAXIT
            END IF
            IF (MAXSIM .NE. MAXSMD) THEN
               WRITE(LUPRI,'(A,I5)')' Maximum vectors in NRCTL     :',
     *                              MAXSIM
            END IF
            IF (NEWRD) THEN
               WRITE (LUPRI,'(/A)') ' New RD file will be used.'
            ELSE IF (NOTRIA) THEN
               WRITE (LUPRI,'(/2A)') ' Old solution vectors not',
     *            ' used as trial vectors in CCZRES.'
            END IF
            IF (NOSOL) THEN
                WRITE (LUPRI,'(/A)') ' Linear equations not solved.'
            END IF
            IF (EXPJAC) THEN
                WRITE (LUPRI,'(/A)') ' Jacobian matrices explicitly'/
     *                /' calculated from left and right are compared.'
            END IF
            IF (LRTEST) THEN
               WRITE (LUPRI,'(/A)') ' Test orthogonality to RHS of '/
     *                /'difference between left and right solutions.'
            END IF
         END IF
      END IF
      RETURN
      END
