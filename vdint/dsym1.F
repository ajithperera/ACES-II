      SUBROUTINE DSYM1(DENMAT,DSCFMAT,FOCMAT,DSO,DSOSCF,
     &                 FSO,NBAST,DENSCF)
C
C     Take density matrix in symmetry orbital basis and generate
C     density matrix over distinct pairs of AOs
C
CEND
C                                          880418  PRT
C  
C  ADAPTED TO THE ACES II ENVIRONEMENT   OCT/90/JG 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL SHARE,DENSCF
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
C
      DIMENSION DSO(1),DSOSCF(1),FSO(1),DENMAT(1),DSCFMAT(1),
     &          FOCMAT(1)
      COMMON /SHELLSi/ KMAX,
     &                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     &                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     &                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     &                NRCO(MXSHEL), NUMCF(MXSHEL),
     &                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3), SHARE(MXSHEL)
      COMMON /PINCOM/ IPIND(MXCORB), IBLOCK(MXCORB), INDGEN(MXCORB)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      COMMON/IPRI/IPRINT
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GT. 10) THEN
       CALL HEADER('Subroutine DSYM1',-1)
      END IF
C
C     Loop over all irreps in molecule
C
      ISOFF = 0
      ISTR = 1
      NNBASX = NBAST*(NBAST + 1)/2
      CALL ZERO(DENMAT,NNBASX)
      CALL ZERO(FOCMAT,NNBASX)
      IF(DENSCF) CALL ZERO(DSCFMAT,NNBASX)
      DO 100 IREP = 0, MAXLOP
         NORBI = NAOS(IREP+1)
         IF (NORBI .EQ. 0) GOTO 110
         DO 200 I = ISTR,ISTR + NORBI - 1
            IA   = IBTSHR(IPIND(I),16)
            NA   = IBTSHR(IPIND(I),8) - IA*256
            IOFF = KSTRT(IA)
            MULA = ISTBAO(IA)
            INDA = IOFF + NA
            DO 300 J = ISTR,I
               IB     = IBTSHR(IPIND(J),16)
               NB     = IBTSHR(IPIND(J),8) - IB*256
               JOFF   = KSTRT(IB)
               NHKTB  = NHKT(IB)
               KHKTB  = KHKT(IB)
               NRCB   = NRCO(IB)
               NBTYP  = (NB-1)/NRCB+1
               MULB   = ISTBAO(IB)
               MAB    = IBTOR(MULA,MULB)
               KAB    = IBTAND(MULA,MULB)
               HKAB   = FMULT(KAB)
               ISOFF  = ISOFF + 1
               DSYMIJ = DSO(ISOFF)
               FSYMIJ = FSO(ISOFF)
               INDB   = JOFF + NB - KHKTB*NRCB
               DO 400 ISYMOP = 0, MAXLOT
                  IF (IBTAND(ISYMOP,MAB) .NE. 0) GOTO 400
                  INDB = INDB + KHKTB*NRCB
C
C                 Weight and parity factor
C
                  FAC = HKAB*
     *                   PT(IBTAND(ISYMOP,IBTXOR(IREP,
     *                       ITYPE(NHKTB,NBTYP))))
                  INDM = MAX(INDA,INDB)
                  IND  = (INDM*(INDM - 3))/2 + INDA + INDB
                  DENMAT(IND) = DENMAT(IND) + FAC*DSYMIJ
                  FOCMAT(IND) = FOCMAT(IND) + FAC*FSYMIJ
400            CONTINUE
               IF(DENSCF) THEN
                DSYMIJ=DSOSCF(ISOFF)
                INDB   = JOFF + NB - KHKTB*NRCB
                DO 410 ISYMOP = 0, MAXLOT
                  IF (IBTAND(ISYMOP,MAB) .NE. 0) GOTO 410
                  INDB = INDB + KHKTB*NRCB
C
C                 Weight and parity factor
C
                  FAC = HKAB*
     *                   PT(IBTAND(ISYMOP,IBTXOR(IREP,
     *                       ITYPE(NHKTB,NBTYP))))
                  INDM=MAX(INDA,INDB)
                  IND=(INDM*(INDM-3))/2+INDA+INDB
                  DSCFMAT(IND)=DSCFMAT(IND)+FAC*DSYMIJ
410             CONTINUE
               ENDIF
300         CONTINUE
200      CONTINUE
110      CONTINUE
         ISTR = ISTR + NORBI
100   CONTINUE
      IPRINT = 0
      IF (IPRINT .GT. 10) THEN
         CALL PRITRI(DENMAT,NBAST,
     *    'Total density matrix (symmetry distinct AO basis)')
         CALL PRITRI(FOCMAT,NBAST,
     *    'Intermediate matrix (symmetry distinct AO basis)')
      END IF
      IPRINT = 0
      RETURN
      END
