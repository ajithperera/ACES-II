      SUBROUTINE PRIHES (SHESS,NCOORD,KEY,ISAVE)
C
C     Adapted from Nelson's OUTPAK for Hessians March 1985 tuh
C     Revised 16-Dec-1983 by Hans Jorgen Aa. Jensen.
C     Jun 28 1988 tuh - modified for symmetry
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (AZERO=0.D00, KCOL=6)
      CHARACTER*(*) KEY
      DIMENSION SHESS(NCOORD,NCOORD), CHESS(MXCOOR,MXCOOR)
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
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      INTEGER BEGIN
      LOGICAL FIRST
      CHARACTER*6 NAMES(MXCOOR), NAME
      SAVE FIRST, NAMES
      DATA FIRST/.TRUE./
C
      IF(ISAVE.EQ.1) THEN
       OPEN(UNIT=82,FILE='FCM',STATUS='UNKNOWN',FORM='FORMATTED')
       REWIND(82)
      ELSEIF(ISAVE.EQ.2) THEN
       OPEN(UNIT=82,FILE='FCM',STATUS='OLD',FORM='FORMATTED')
       REWIND(82)
      ENDIF
      IF (FIRST) THEN
         FIRST = .FALSE.
         IOFF = 0
         JOFF = 0
         DO 10 IATOM = 1, NUCIND
            MULTI = MULT(ISTBNU(IATOM))
            IF (MULTI .EQ. 1) THEN
               NAMES(IOFF + 1) = NAMEX(JOFF + 1)
               NAMES(IOFF + 2) = NAMEX(JOFF + 2)
               NAMES(IOFF + 3) = NAMEX(JOFF + 3)
               IOFF = IOFF + 3
            ELSE
               DO 15 J = 1, MULTI
                  NAMES(IOFF + 1) = NAMEX(JOFF + 1)(1:3)//CHAR(48+J)
     *                            //NAMEX(JOFF + 1)(5:6)
                  NAMES(IOFF + 2) = NAMEX(JOFF + 2)(1:3)//CHAR(48+J)
     *                            //NAMEX(JOFF + 2)(5:6)
                  NAMES(IOFF + 3) = NAMEX(JOFF + 3)(1:3)//CHAR(48+J)
     *                            //NAMEX(JOFF + 3)(5:6)
                  IOFF = IOFF + 3
   15          CONTINUE
            END IF
            JOFF = JOFF + 3
   10    CONTINUE
      END IF
C
      IF (KEY .EQ. 'CENTERS') THEN
         NROW = 3*NUCDEP
      ELSE IF (KEY .EQ. 'ATOMS') THEN
         NROW = 3*NATOMS
      ELSE
         WRITE (LUPRI,'(//A/2A//)') ' >>> ERROR in PRIHES ',
     *      '     - incorrect specification of keyword:', KEY
         STOP 
      END IF
      IF (MAXLOP .EQ. 0) THEN
         LAST = MIN(NROW,KCOL)
         BEGIN= 1
   50    CONTINUE
cmn         WRITE (LUPRI,1000) (NAMEX(I),I = BEGIN,LAST)
cmn         WRITE (LUPRI,'()')
         NCOL = 1
         DO 100 K = BEGIN,NROW
            DO 200 I = 1,NCOL
               IF (SHESS(K,(BEGIN-1)+I) .NE. AZERO) GO TO 400
  200       CONTINUE
            GO TO 300
 400        CONTINUE
cmn         WRITE (LUPRI,2000)' ',NAMEX(K),
cmn     *            (SHESS(K,(BEGIN-1)+J),J=1,NCOL)
cmn            IF (MOD(K,3) .EQ. 0) WRITE (LUPRI,'()')
  300       IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
  100    CONTINUE
cmn         WRITE (LUPRI,'()')
         LAST = MIN(LAST+KCOL,NROW)
         BEGIN= BEGIN+NCOL
         IF (BEGIN.LE.NROW) GO TO 50
      ELSE
         NOFF = 0
         DO 500 IREP = 0, MAXLOP
         IF (NCRREP(IREP) .GT. 0) THEN
cmn            WRITE (LUPRI,'(/11X,A,I1/)') 'Symmetry ',IREP+1
            NROW  = NOFF + NCRREP(IREP)
            LAST  = MIN(NROW,NOFF + KCOL)
            BEGIN = NOFF + 1
C
  550       CONTINUE
cmn            WRITE (LUPRI,1000) (NAMEX(IPTCOR(I)),I = BEGIN,LAST)
cmn            WRITE (LUPRI,'()')
            NCOL = 1
            DO 40 K = BEGIN,NROW
 20            continue
cmn   20          WRITE (LUPRI,2000) ' ', NAMEX(IPTCOR(K)),
cmn     *               (SHESS(K,(BEGIN-1)+J),J=1,NCOL)
   30          IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
   40       CONTINUE
cmn            WRITE (LUPRI,'()')
            LAST = MIN(LAST+KCOL,NROW)
            BEGIN= BEGIN+NCOL
            IF (BEGIN.LE.NROW) GO TO 550
C
            NOFF = NOFF + NCRREP(IREP)
         END IF
  500    CONTINUE
C
C        Print in non-symmetry basis
C
         IF (.TRUE.) THEN
cmn            WRITE (LUPRI,'(//)')
            CALL ZERO(CHESS,MXCOOR*MXCOOR)
            CALL TRAHES(SHESS,CHESS,NCOORD)
            NROW = NCOORD
            LAST = MIN(NROW,KCOL)
            BEGIN= 1
   51       CONTINUE
cmn            WRITE (LUPRI,1000) (NAMES(I),I = BEGIN,LAST)
cmn            WRITE (LUPRI,'()')
            NCOL = 1
            DO 101 K = BEGIN,NROW
               DO 201 I = 1,NCOL
                  IF (CHESS(K,(BEGIN-1)+I) .NE. AZERO) GO TO 401
  201          CONTINUE
               GO TO 301
 401            continue
cmn  401          WRITE (LUPRI,2000)' ',NAMES(K),
cmn     *               (CHESS(K,(BEGIN-1)+J),J=1,NCOL)
cmn               IF (MOD(K,3) .EQ. 0) WRITE (LUPRI,'()')
  301          IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
  101       CONTINUE
cmn            WRITE (LUPRI,'()')
            LAST = MIN(LAST+KCOL,NROW)
            BEGIN= BEGIN+NCOL
            IF (BEGIN.LE.NROW) GO TO 51
         END IF
      END IF
cmn      WRITE (LUPRI,'()')
      IF(ISAVE.NE.0) THEN
      IF(MAXLOP.EQ.0) THEN
       IF(ISAVE.EQ.2) THEN
        READ(82,*) IJUNK,JJUNK
        DO 590 I=1,NCOORD
         READ(82,3001)(CHESS(I,J),J=1,NCOORD)
590     CONTINUE
        REWIND(82)
        DO 595 I=1,NCOORD
        DO 595 J=1,NCOORD
         SHESS(I,J)=SHESS(I,J)+CHESS(I,J)
595     CONTINUE
       ENDIF
       WRITE(82,3000) NCOORD/3,NCOORD*2
3000   FORMAT(2I5)
       DO 600 I=1,NCOORD
        WRITE(82,3001)(SHESS(I,J), J=1,NCOORD)
3001    FORMAT(3F20.10)
600    CONTINUE
       CALL PUTREC(20,'JOBARC','HESSIANM',NCOORD*NCOORD*IINTFP,SHESS)
      ELSE
       IF(ISAVE.EQ.2) THEN
        READ(82,*) IJUNK,JJUNK
        DO 690 I=1,NCOORD
         READ(82,3001)(SHESS(I,J),J=1,NCOORD)
690     CONTINUE
        REWIND(82)
        DO 695 I=1,NCOORD
        DO 695 J=1,NCOORD
         CHESS(I,J)=SHESS(I,J)+CHESS(I,J)
695     CONTINUE
       ENDIF
       WRITE(82,3000) NCOORD/3,NCOORD*2
       DO 700 I=1,NCOORD
       WRITE(82,3001)(CHESS(I,J), J=1,NCOORD)
700    CONTINUE
       CALL PUTREC(20,'JOBARC','HESSIANM',NCOORD*NCOORD*IINTFP,CHESS)
      ENDIF
      CLOSE(UNIT=82,STATUS='KEEP')
      ENDIF
      RETURN
 1000 FORMAT (8X,6(3X,A6,3X),(3X,A6,3X))
 2000 FORMAT (A1,A6,6F12.6)
      END
