      SUBROUTINE DSYM1A(DEN,DSOA,DSOB,NBAST,IPRINT)
C
C     Take density matrix in symmetry orbital basis and generate
C     density matrix over distinct pairs of AOs
C
CEND
C                                          880418  PRT
C  
C  ADAPTED TO HE ACES II ENVIRONEMENT   OCT/90/JG 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      DIMENSION DSOA(1), DSOB(1),DEN(1)
      LOGICAL SHARE
      LOGICAL SCF,NONHF
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
      COMMON/METHOD/IUHF,SCF,NONHF
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine DSYM1',-1)
      END IF
C
C     Loop over all irreps in molecule
C
      ISOFF = 0
      ISTR = 1
      NNBASX = NBAST*(NBAST + 1)/2
      CALL ZERO(DSOA,NNBASX)
      IOFFB=NBASTT*IUHF
      IF(IUHF.NE.0) CALL ZERO(DSOB,NNBASX)
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
               DSYMIJA = DEN(ISOFF)
               DSYMIJB = DEN(ISOFF+IOFFB)
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
                  DSOA(IND) = DSOA(IND) + FAC*DSYMIJA
                  IF(IUHF.NE.0) DSOB(IND) = DSOB(IND) + FAC*DSYMIJB
400            CONTINUE
300         CONTINUE
200      CONTINUE
110      CONTINUE
         ISTR = ISTR + NORBI
100   CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL PRITRI(DSOA,NBAST,
     *    'Folded alpha density matrix (symmetry distinct AO basis)')
        IF(IUHF.NE.0) THEN
         CALL PRITRI(DSOB,NBAST,
     *    'Folded beta density matrix (symmetry distinct AO basis)')
       ENDIF
      ENDIF
      RETURN
      END
