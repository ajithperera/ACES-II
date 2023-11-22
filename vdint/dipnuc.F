      SUBROUTINE DIPNUC(DIFDIP,DIPMN,DDIPN,NCOORD)
C
C     CALCULATES NUCLEAR CONTRIBUTIONS TO ELECTRIC DIPOLE MOMENTS
C     AND DIPOLE GRADIENTS
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (D0 = 0.0D0, D100 = 100.0D0)
      LOGICAL DIFDIP,ELECT,MAGNET,GEOM,MAGNET2,GIAO
      LOGICAL MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *        VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *        H2MO
C
      DIMENSION DIPMN(3),DDIPN(3,NCOORD)
      COMMON /ABAINF/ IPRDEF,
     *                MOLGRD, MOLHES, DIPDER, POLAR,  INPTES,
     *                VIB,    RESTAR, DOWALK, GDALL,  FOCK,
     *                H2MO
      CHARACTER NAMEX*6
      LOGICAL DCORD, DCORGD, NOORBT, DOPERT
      COMMON/PROP/ELECT,MAGNET,GEOM,MAGNET2,GIAO
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
      COMMON/IPRI/IPRINT
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     Dipole moment
C
      CALL ZERO(DIPMN,3)
      DO 100 IATOM = 1, NUCIND
         CHA = CHARGE(IATOM)
         IF (ABS(CHA) .LT. D100) THEN
            FAC = FMULT(ISTBNU(IATOM))*CHA
            DO 200 ICOOR = 1, 3
               IF (ISYTYP(ICOOR) .EQ. 0) THEN
                  DIPMN(ICOOR) = DIPMN(ICOOR) + FAC*CORD(IATOM,ICOOR)
               END IF
  200       CONTINUE
         END IF
  100 CONTINUE
C
C     Dipole moment gradient
C
      IF (DIFDIP.AND.(.NOT.ELECT)) THEN
         CALL ZERO(DDIPN,3*NCOORD)
         DO 300 IREP = 0, MAXLOP
            IF (NAXREP(IREP) .GT. 0) THEN
               DO 400 IATOM = 1, NUCIND
                  CHA = CHARGE(IATOM)
                  IF (ABS(CHA) .LT. D100) THEN
                     FAC = FMULT(ISTBNU(IATOM))*CHA
                     DO 500 ICOOR = 1, 3
                        IF (ISYTYP(ICOOR) .EQ. IREP) THEN
                           ISC = IPTCNT(3*(IATOM - 1) + ICOOR,IREP)
                           IF (ISC.GT.0) DDIPN(IPTAX(ICOOR),ISC) = FAC
                        END IF
  500                CONTINUE
                  END IF
  400          CONTINUE
            END IF
  300    CONTINUE
      END IF
      IF (IPRINT .GT. 0) THEN
       CALL HEADER('Nuclear contribution to dipole moments',-1)
       CALL DP0PRI(DIPMN)
       IF (DIFDIP.AND.(.NOT.ELECT)) THEN
        CALL HEADER('Nuclear contribution to dipole derivatives',-1)
        CALL DP1PRI(DDIPN,NCOORD,'CENTERS',0)
       END IF
      END IF
      RETURN
      END
