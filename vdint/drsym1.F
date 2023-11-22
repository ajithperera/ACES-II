      SUBROUTINE DRSYM1(SAODER,TAODER,WORK1,KSDAO,KFDAO,
     *                  ICENTA,ICENTB,ISYMOP,MULA,MULB,NHKTA,NHKTB,
     *                  KHKTA,KHKTB,HKAB,LDIAG,LAEQB,
     *                  THRESH,IPRINT)
C
C     Arrange calculation of symmetry-adapted integral derivatives
C     from distinct AO integral derivatives
C                                                PRT & TUH  880428
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      DIMENSION SAODER(1), TAODER(1),WORK1(1)
      DIMENSION KSDAO(8),KFDAO(8)
      LOGICAL LDIAG, FULMAT, LAEQB
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
      IBTXOR(I,J) = XOR(I,J)
      IF (IPRINT .GT. 10) THEN
       CALL HEADER('Subroutine DRSYM1',-1)
       WRITE (LUPRI,'(A,2I5)') ' ICENTA/B ', ICENTA, ICENTB
       WRITE (LUPRI,'(A, I5)') ' ISYMOP   ', ISYMOP
       WRITE (LUPRI,'(A,2I5)') ' NHKTA/B  ', NHKTA, NHKTB
       WRITE (LUPRI,'(A,2I5)') ' KHKTA/B  ', KHKTA, KHKTB
       WRITE (LUPRI,'(A,2I5)') ' MULA/B   ', MULA, MULB
       WRITE (LUPRI,'(A,F12.6)') ' HKAB   ', HKAB
       WRITE (LUPRI,'(A, L5)') ' LDIAG    ', LDIAG
      END IF
      FULMAT = .TRUE.
      NMATS = 3*NUCDEP*(MAXLOP+1)
      KHKTAB = KHKTA*KHKTB
      NRCAB=NRCA*NRCB
      IF(LDIAG) NRCAB=NRCA*(NRCA+1)/2
      NINT=KHKTAB*NRCAB
C
C LOOP OVER RIGHT AND LEFT HAND SIDE OF THE 1E INTEGRAL
C
      DO 10 ICL = 1,2
C
C  Determine factors to account for use of transl. invariance
C
       IF (ICL .EQ. 1) THEN
        JCENT = ICENTA
       ELSE
        JCENT = ICENTB
       ENDIF
       MULJ = ISTBNU(JCENT)
       IOFF = 1
C
C  Run over Cartesian directions
C
       DO 20 JCDIR = 1,3
        IDIR = JCDIR - 1
        ISYTYJ = ISYTYP(JCDIR)
C
C Run over irrep's of the differentiation operator
C
        DO 30 IREPD = 0,MAXLOP
         IF (IBTAND(MULJ,IBTXOR(IREPD,ISYTYJ)) .EQ. 0) THEN
          IF (ICL .EQ. 1) THEN
           FAC = HKAB
          ELSE
           FAC = - HKAB*PT(IBTAND(ISYTYJ,ISYMOP))
     *           *PT(IBTAND(IREPD, ISYMOP))
          END IF
          IF (IPRINT .GT. 20) THEN
           WRITE (LUPRI, '(/A,3I5)') ' ICL, JCDIR, IREPD ',
     *                                 ICL, JCDIR, IREPD
           WRITE (LUPRI, '(A,2I5)') ' MULJ, ISYTYJ ', MULJ, ISYTYJ
           WRITE (LUPRI, '(A,F12.6)') ' FAC ', FAC
          END IF
          IMAT0 = (IPTCNT(3*(JCENT-1)+IDIR+1,IREPD)-1)*(MAXLOP+1)+1
          IF (IREPD .EQ. 0) THEN
C
C  Overlap matrix - totally symmetric perturbation
C
c         IF(ICL.EQ.1) THEN
          IOFFSET=(IPTCNT(3*(JCENT-1)+IDIR+1,0)-1)*NLENT(1)
          CALL SYM1S(SAODER(IOFF),WORK1(KSDAO(1)+IOFFSET),
     *               ISYMOP,MULA,MULB,
     *               NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,LAEQB,
     *               FULMAT,
     *               THRESH,IMAT0,1,IPRINT,ONE)
c         ENDIF
C
C  Kinetic energy  - totally symmetric perturbation
C
          IMAT0 = NMATS + IMAT0
          IOFFSET=(IPTCNT(3*(JCENT-1)+IDIR+1,0)-1)*NLENT(1)
          CALL SYM1S(TAODER(IOFF),WORK1(KFDAO(1)+IOFFSET),
     *               ISYMOP,MULA,MULB,
     *               NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,LAEQB,
     *               FULMAT,
     *               THRESH,IMAT0,1,IPRINT,ONE)
         ELSE
C
C  Overlap matrix - non-symmetric perturbation
C
c              IF(ICL.EQ.1) THEN
          IOFFSET=(IPTCNT(3*(JCENT-1)+IDIR+1,IREPD)
     *            -IPERT(IREPD+1)-1)
     *            *NLENT(IREPD+1)
          CALL SYM1N(SAODER(IOFF),
     *               WORK1(KSDAO(IREPD+1)+IOFFSET),
     *               IREPD,ISYMOP,MULA,
     *               MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *               LAEQB,FULMAT,THRESH,IMAT0,1,IPRINT,ONE)
c        ENDIF
C
C  Kinetic energy  - non-symmetric perturbation
C
c        IMAT0 = NMATS + IMAT0
         IOFFSET=(IPTCNT(3*(JCENT-1)+IDIR+1,IREPD)
     *           -IPERT(IREPD+1)-1)
     *           *NLENT(IREPD+1)
         CALL SYM1N(TAODER(IOFF),
     *              WORK1(KFDAO(IREPD+1)+IOFFSET),
     *              IREPD,ISYMOP,MULA,
     *              MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *              LAEQB,FULMAT,THRESH,IMAT0,1,IPRINT,ONE)
        END IF
       END IF
30    CONTINUE
      IOFF = IOFF + NINT
20    CONTINUE
10    CONTINUE
      RETURN
      END
