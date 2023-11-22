      SUBROUTINE DSYM2A(DHFSO,NBAST,IPRINT,DSOA,DSOB,NIR)
C
C     Take total and active density matrices in symmetry
C     orbital basis and generate density matrices (non-folded)
C     over distinct pairs of AOs
C
CEND
C                                          880418  PRT
C
C  ADAPTED TO THE ACES II ENVIRONMENT AND EXTENDED FOR UHF AND NON-HF
C  OCT/90/JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL SHARE
c&line mod
      LOGICAL SCF,NONHF
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      DIMENSION DHFSO(1)
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
      DIMENSION DSOA(NBAST*(NBAST+1)/2,0:NIR),
     *          DSOB(NBAST*(NBAST+1)/2,0:NIR)
      COMMON/METHOD/IUHF,SCF,NONHF
c&line del
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
C
      DATA ONE,TWO /1.0D0,2.D0/
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine DSYM2A',-1)
      END IF
C
      IF(IUHF.EQ.0) THEN
       CALL ZERO(DSOA(1,0),NBAST*(NBAST+1)/2*(NIR+1))
      ELSE
       CALL ZERO(DSOA(1,0),NBAST*(NBAST+1)/2*(NIR+1))
       CALL ZERO(DSOB(1,0),NBAST*(NBAST+1)/2*(NIR+1))
      ENDIF
C
      IOFFA=0
      IOFFB=IOFFA+NBASTT*IUHF
C
C     Loop over all irreps in molecule
C
      ISOFF = 0
      ISTR = 1
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
               IB    = IBTSHR(IPIND(J),16)
               NB    = IBTSHR(IPIND(J),8) - IB*256
               JOFF  = KSTRT(IB)
               NHKTB = NHKT(IB)
               KHKTB = KHKT(IB)
               MULB  = ISTBAO(IB)
               MAB   = IBTOR(MULA,MULB)
               ISOFF = ISOFF + 1
               DHFAIJ = DHFSO(ISOFF)
               DHFBIJ = DHFSO(ISOFF+IOFFB)
               DO 400 ISYMOP = 0, MAXLOT
                  FAC = PT(IBTAND(ISYMOP,IREP))
                  INDB = JOFF + NB
                  INDM=MAX(INDA,INDB)
                  IND=(INDM*(INDM-3))/2+INDA+INDB
                  DSOA(IND,ISYMOP) = DSOA(IND,ISYMOP)
     *                                     + FAC*DHFAIJ
                IF(IUHF.NE.0) THEN
                  DSOB(IND,ISYMOP) = DSOB(IND,ISYMOP)
     *                                     + FAC*DHFBIJ
                ENDIF
400            CONTINUE
300         CONTINUE
200      CONTINUE
110      CONTINUE
         ISTR = ISTR + NORBI
100   CONTINUE
      IPRINT = 0
      IF (IPRINT .GT. 10) THEN
         DO 500 ISYMOP = 0, MAXLOP
          CALL PRITRI(DSOA(1,ISYMOP),NBAST,
     *        'Folded alpha density matrix (symmety distint AO basis)')
          IF(IUHF.NE.0) THEN
          CALL PRITRI(DSOB(1,ISYMOP),NBAST,
     *        'Folded alpha density matrix (symmety distint AO basis)')
          ENDIF
  500    CONTINUE
      END IF
      IPRINT = 0
C
      RETURN
      END
