      SUBROUTINE DASYM4(CAODER,WORK1,KFDAO,NATOMC,ICENTC,
     &                  ISYMOP,JSYMOP,MULA,MULB,NHKTA,NHKTB,
     &                  KHKTA,KHKTB,HKAB,LDIAG,LAEQB,
     &                  THRESH,IPRINT)
C
C ARRANGE CALCULATION OF SYMMETRY-ADAPTED INTEGRAL DERIVATIVES
C FROM DISTINCT AO SPIN DIPOLAR INTEGRALS
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      LOGICAL LDIAG, ONECEN, LAEQB
      DIMENSION CAODER(1000),WORK1(1000),KFDAO(8)
      LOGICAL FULMAT
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
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
C
      COMMON/GENCON/NRCA,NRCB,CONTA(MXCONT),CONTB(MXCONT)
      COMMON/LSYM/NLENQ(8),NLENT(8)
      COMMON/PERT/NTPERT,NPERT(8),IPERT(8),IXPERT,IYPERT,IZPERT,
     &            IYZPERT,IXZPERT,IXYPERT,ITRANSX,ITRANSY,ITRANSZ,
     &            NUCIND1
C
      DATA ONE /1.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      FULMAT = .TRUE.
      KHKTAB = KHKTA*KHKTB
      NRCAB=NRCA*NRCB
      IF(LDIAG) NRCAB=NRCA*(NRCA+1)/2
      NINTS=KHKTAB*NRCAB*NATOMC
      NMATS = 3*NUCDEP*(MAXLOP+1)
C
      JCENT=ICENTC
      MULJ = ISTBNU(JCENT)
      IOFF = 1
C
C RUN OVER CARTESIAN DIRECTIONS
C
      IJDIR=0
      DO 20 JCDIR =1,3
       IDIR=JCDIR-1
       DO 20 KCDIR=JCDIR,3
        JDIR=JCDIR-1
         IJDIR=IJDIR+1
         ISYTIJ=IBTXOR(ISYTYP(JCDIR),ISYTYP(KCDIR))
C
C RUN OVER IRREPS OF THE DIFFERENTIATION OPERATOR
C
         DO 30 IREPD = 0,MAXLOP
          IF (IBTAND(MULJ,IBTXOR(IREPD,ISYTIJ)) .EQ. 0) THEN
           FAC = HKAB*PT(IBTAND(ISYTIJ,JSYMOP))
     &               *PT(IBTAND(IREPD ,JSYMOP))
           IOFFSET=(IPTCNT(6*(JCENT-1)+IJDIR,IREPD)
     &             -IPERT(IREPD+1)-1)*NLENT(IREPD+1)
           IF (IREPD .EQ. 0) THEN
            CALL SYM1S(CAODER(IOFF),WORK1(KFDAO(IREPD+1)
     &                 +IOFFSET),ISYMOP,MULA,MULB,
     &                 NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,LAEQB,
     &                 FULMAT,
     &                 THRESH,IMAT0,NATOMC,IPRINT,ONE)
           ELSE
            CALL SYM1N(CAODER(IOFF),WORK1(KFDAO(IREPD+1)
     &                 +IOFFSET),IREPD,ISYMOP,MULA,
     &                 MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     &                 LAEQB,
     &                 FULMAT,THRESH,IMAT0,NATOMC,IPRINT,ONE)
           END IF
          END IF
30       CONTINUE
         IOFF = IOFF + NINTS
20     CONTINUE
10    CONTINUE
C
      RETURN
      END
