      SUBROUTINE CSHPRI(AMAT,NCOORD,ISAVE,TOTAL)
C
C PRINT ROUTINE FOR SHIELDING TENSORS
C
CEND
C
C SEPTEMBER/91   JG  UNIVERSITY OF KARLSRUHE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      CHARACTER*1 CHRXYZ(-3:3)
      CHARACTER*7 KEY
      DIMENSION AMAT(3,NCOORD), CMAT(3,MXCOOR),MULTAT(MXCENT)
      CHARACTER NAMEX*6
      LOGICAL DCORD,DCORGD,NOORBT,DOPERT,TOTAL
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
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA CHRXYZ /'z','y','x',' ','X','Y','Z'/
      DATA AZERO,THREE /0.D0,3.D0/
C
      IF(TOTAL) THEN
       IF(ISAVE.EQ.1) THEN
        OPEN(UNIT=82,FILE='CSHIFT',STATUS='UNKNOWN',FORM='FORMATTED')
        REWIND(82)
       ELSEIF(ISAVE.EQ.2) THEN
        OPEN(UNIT=82,FILE='CSHIFT',STATUS='OLD',FORM='FORMATTED')
        REWIND(82)
       ENDIF
      ELSE
       IF(ISAVE.EQ.1) THEN
        OPEN(UNIT=82,FILE='CSHIFTSCF',STATUS='UNKNOWN',
     &       FORM='FORMATTED')
        REWIND(82)
       ELSEIF(ISAVE.EQ.2) THEN
        OPEN(UNIT=82,FILE='CSHIFTSCF',STATUS='OLD',
     &       FORM='FORMATTED')
        REWIND(82)
       ENDIF
      ENDIF
      DO 1 INUC=1,NUCIND
       MULTAT(INUC)=1
1      CONTINUE
      IOFFAX = 0
      DO 100 IREP = 0, MAXLOP
       NAXIS = NAXREP(IREP)
       IF ((NAXIS .GT. 0) .AND. (NCRREP(IREP) .GT. 0)) THEN
        IF (MAXLOP .GT. 0) THEN
         WRITE (LUPRI,'(/11X,A,I1)') 'Symmetry ', IREP + 1
        END IF
        WRITE (LUPRI,'(30X,3(2A,13X))')
     *  ('B',CHRXYZ(-IPTXYZ(I,IREP)),I=1,NAXIS)
        WRITE (LUPRI,'()')
        DO 200 IATOM = 1, NUCIND
c         AVER=AZERO
         DO 300 ICOOR = 1, 3
          ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP)
          IF (ISCOOR .GT. 0) THEN
           WRITE (LUPRI,'(15X,A6,3F15.8)')
     *         NAMEX(IPTCOR(ISCOOR)),
     *         (AMAT(IOFFAX + I,ISCOOR), I=1,NAXIS)
          END IF
  300    CONTINUE
         IF (MAXLOP .EQ. 0) THEN
          WRITE (LUPRI,'()')
         ENDIF
  200   CONTINUE
       END IF
       IOFFAX = IOFFAX + NAXIS
  100 CONTINUE
      WRITE (LUPRI,'(/)')
      IF(MAXLOP.EQ.0) THEN

       IF(ISAVE.NE.0) THEN
        IF(ISAVE.EQ.2) THEN
         DO 319 IATOM=1,NUCIND
          DO 309 IFIELD=1,3
           READ(82,5101)CMAT(IFIELD,3*IATOM-2),
     &          CMAT(IFIELD,3*IATOM-1),CMAT(IFIELD,3*IATOM)
           AMAT(IFIELD,3*IATOM-2)=CMAT(IFIELD,3*IATOM-2)
     &                           +AMAT(IFIELD,3*IATOM-2)
           AMAT(IFIELD,3*IATOM-1)=CMAT(IFIELD,3*IATOM-1)
     &                           +AMAT(IFIELD,3*IATOM-1)
           AMAT(IFIELD,3*IATOM)=CMAT(IFIELD,3*IATOM)
     &                           +AMAT(IFIELD,3*IATOM)
309       CONTINUE
319      CONTINUE
         REWIND(82)
        ENDIF
        DO 320 IATOM=1,NUCIND
         DO 310 IFIELD=1,3
          WRITE(82,5101)AMAT(IFIELD,3*IATOM-2),
     &         AMAT(IFIELD,3*IATOM-1),AMAT(IFIELD,3*IATOM)
5101      FORMAT(4F20.10)
310      CONTINUE
320     CONTINUE
       ENDIF
      ELSE
       NCOOR = 3*NUCDEP
       CALL TRACHS(AMAT,CMAT,NCOOR)
       WRITE (LUPRI,'()')
       WRITE (LUPRI,'(13X,3(10X,A2))') 'Bx', 'By', 'Bz'
       WRITE (LUPRI,'()')
       IOFF = 1
       DO 400 IATOM = 1, NUCIND
        MULTI = MULT(ISTBNU(IATOM))
        MULTAT(IATOM)=MULTI
        ISTR  = 3*IATOM - 3
        IF (MULTI .EQ. 1) THEN
         DO 500 ICOOR = 1, 3
          WRITE (LUPRI,1000) ' ', NAMEX(ISTR+ICOOR)(1:4),' ',
     *           NAMEX(ISTR+ICOOR)(6:6), (CMAT(K,IOFF),K=1,3)
          IOFF = IOFF + 1
  500    CONTINUE
         WRITE (LUPRI,'()')
        ELSE
         DO 600 J = 1, MULTI
          DO 700 ICOOR = 1, 3
           WRITE (LUPRI,1010) ' ',NAMEX(ISTR+ICOOR)(1:4),J,
     *            NAMEX(ISTR+ICOOR)(6:6), (CMAT(K,IOFF),K=1,3)
                  IOFF = IOFF + 1
  700     CONTINUE
          WRITE (LUPRI,'()')
  600    CONTINUE
        END IF
  400  CONTINUE
       IF(ISAVE.NE.0) THEN
        IF(ISAVE.EQ.2) THEN
         DO 419 IATOM=1,NUCIND
          DO 409 IFIELD=1,3
           READ(82,5101)AMAT(IFIELD,3*IATOM-2),
     &         AMAT(IFIELD,3*IATOM-1),AMAT(IFIELD,3*IATOM)
               CMAT(IFIELD,3*IATOM-2)=CMAT(IFIELD,3*IATOM-2)
     &                           +AMAT(IFIELD,3*IATOM-2)
             CMAT(IFIELD,3*IATOM-1)=CMAT(IFIELD,3*IATOM-1)
     &                           +AMAT(IFIELD,3*IATOM-1)
             CMAT(IFIELD,3*IATOM)=CMAT(IFIELD,3*IATOM)
     &                           +AMAT(IFIELD,3*IATOM)
409         CONTINUE
419        CONTINUE
           REWIND(82)
          ENDIF
          JATOM=0
          DO 420 IATOM=1,NUCIND
           DO 420 J=1,MULT(ISTBNU(IATOM))
           JATOM=JATOM+1
           DO 410 IFIELD=1,3
            WRITE(82,5101)CMAT(IFIELD,3*JATOM-2),
     &           CMAT(IFIELD,3*JATOM-1),CMAT(IFIELD,3*JATOM)
410      CONTINUE
420     CONTINUE
       ENDIF
      END IF
      IF(ISAVE.NE.0) CLOSE(UNIT=82,STATUS='KEEP')
      IF(TOTAL) THEN
       call putrec(20,'JOBARC','NINDATOM',1,NUCIND)
       call putrec(20,'JOBARC','MULTATOM',NUCIND,MULTAT)
       call putrec(20,'JOBARC','NAMCOORD',IINTFP*3*NUCIND,NAMEX)
      ENDIF 
 1000 FORMAT (5X,A1,A4,A1,3X,A1,3F12.6)
 1010 FORMAT (5X,A1,A4,I1,3X,A1,3F12.6)
      RETURN
      END
