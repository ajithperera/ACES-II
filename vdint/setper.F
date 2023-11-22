      SUBROUTINE SETPER(IPRINT)
C
C     Sets DCORD, DCORGD, NGDVEC, IGDVEC, IPTGDV, NDCOOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
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
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
C
C     DCORD & DCORGD
C
      DOPERT(0) = .FALSE.
      DO 100 JATOM = 1, NUCIND
         DO 110 ICOOR = 1, 3
            JCOOR = 3*(JATOM - 1) + ICOOR
            DCORD (JATOM,ICOOR) = .FALSE.
            DCORGD(JATOM,ICOOR) = .FALSE.
            DO 120 IREP = 0, MAXLOP
               ISYMCR = IPTCNT(JCOOR,IREP)
               IF (ISYMCR.NE.0) THEN
                  DCORD (JATOM,ICOOR) = DCORD (JATOM,ICOOR)
     *                                  .OR. DOPERT(ISYMCR)
                  DCORGD(JATOM,ICOOR) = DCORGD(JATOM,ICOOR)
     *                                  .OR. DOPERT(ISYMCR)
               END IF
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 30) THEN
         CALL HEADER('DCORD',-1)
         DO 200 JATOM = 1, NUCIND
            WRITE (LUPRI, '(1X,3L5)') (DCORD(JATOM,I),I=1,3)
  200    CONTINUE
         CALL HEADER('DCORGD',-1)
         DO 210 JATOM = 1, NUCIND
            WRITE (LUPRI, '(1X,3L5)') (DCORGD(JATOM,I),I=1,3)
  210    CONTINUE
      END IF
C
C     Number of coordinates in each symmetry for which linear equations
C     will be solved
C
C     NGDVEC(8)        - # GD vectors in this symmetry
C     IGDVEC(8)        - # GD of vectors in previous symmetries
C     IPTGDV(3*NUCIND) - points from GD vector to symmetry coordinate
C
      NOFF = 0
      DO 300 IREP = 0, MAXLOP
         NSYMCR = 0
         IGDVEC(IREP + 1) = NOFF
         DO 310 JATOM = 1, NUCIND
            DO 320 ICOOR = 1, 3
               JCOOR = 3*(JATOM - 1) + ICOOR
               ISCOOR = IPTCNT(JCOOR,IREP)
               IF (ISCOOR .NE. 0) THEN
                  IF (DOPERT(ISCOOR)) THEN
                     NOFF = NOFF + 1
                     NSYMCR = NSYMCR + 1
                     IPTGDV(NOFF) = IPTCNT(JCOOR,IREP)
                  END IF
               END IF
 320        CONTINUE
 310     CONTINUE
         NGDVEC(IREP + 1) = NSYMCR
 300  CONTINUE
      IF (IPRINT .GT. 30) THEN
         WRITE (LUPRI,'(//1X,A,8I5)')' NGDVEC ',(NGDVEC(I),I=1,MAXLOP+1)
         WRITE (LUPRI,'(1X,A,8I5)')  ' IGDVEC ',(IGDVEC(I),I=1,MAXLOP+1)
         WRITE (LUPRI,'(1X,A,9I5,(/9I5))')
     *      ' IPTGDV ', (IPTGDV(I),I=1,NOFF)
      END IF
C
      NDCORD = 0
      DO 400 JATOM = 1,NUCIND
         DO 410 ICOOR = 1, 3
            IF (DCORD(JATOM,ICOOR)) THEN
               NDCORD = NDCORD + 1
               NDCOOR(NDCORD) = 3*(JATOM - 1) + ICOOR
            END IF
  410    CONTINUE
  400 CONTINUE
      RETURN
      END
