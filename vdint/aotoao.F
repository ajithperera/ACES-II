      SUBROUTINE AOTOAO(IAOAO,IPRINT)
C
C     This subroutine sets up pointer IAOAO which convertes between
C     two different orderings of AO's: from the ordering
C     of MOLECULE to an ordering in which the outer loop is over
C     atoms. (This orderings are identical when no symmetry is used.)
C
C     The purpose of this ordering is to make the AO's appear in the
C     same order as in the corresponding run with no symmetry.
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
C
      DIMENSION IAOAO(MXCORB)
C
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
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
      LOGICAL SHARE
      COMMON /SHELLSi/ KMAX,
     &                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     &                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     &                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     &                NRCO(MXSHEL), NUMCF(MXSHEL),
     &                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3), SHARE(MXSHEL)
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GE. 10) CALL TITLER('Output from AOTOAO','*',103)
C
C     Loop over all atoms
C
      IAO   = 0
      JATOM = 0
      DO 100 IATOM = 1, NUCIND
         ISTABA = ISTBNU(IATOM)
         DO 200 ISYMOP = 0, MAXLOP
         IF (IBTAND(ISYMOP,ISTABA) .EQ. 0) THEN
            JATOM = JATOM + 1
C
C           Loop over all orbitals
C
            JAO = 0
            DO 300 ISHELL = 1, KMAX
               ICENT  = NCENT(ISHELL)
               ISTABO = ISTBAO(ISHELL)
               DO 400 ICMP = 1, KHKT(ISHELL)
                do 400 irc=1,nrco(ishell)
                  DO 500 JSYMOP = 0, MAXLOP
                  IF (IBTAND(JSYMOP,ISTABO) .EQ. 0) THEN
                     JAO = JAO + 1
                     JCENT = NUCNUM(ICENT,JSYMOP+1)
                     IF (JATOM .EQ. JCENT) THEN
                        IAO = IAO + 1
                        IAOAO(JAO) = IAO
                     END IF
                  END IF
  500             CONTINUE
  400          CONTINUE
  300       CONTINUE
C
C           End loop over orbitals
C
         END IF
  200    CONTINUE
  100 CONTINUE
C
C     End loop over atoms
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('I - IAOAO',6)
         DO 600 I = 1, IAO
            WRITE (LUPRI,'(4X,I5,1X,I5)') I, IAOAO(I)
  600    CONTINUE
      END IF
      RETURN
      END
