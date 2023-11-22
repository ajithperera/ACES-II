      SUBROUTINE SYM1S(AO,SO,KB,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     &                 HKAB,LDIAG,LAEQB,FULMAT,THRESH,IMAT0,
     &                 JSKIP,IPRINT,ASYM)
C
C     TAKE BLOCK OF DISTINCT AO INTEGRAL (DERIVATIVES) AND 
C     GENERATE SYMMETRIZED CONTRIBUTIONS TO SO INTEGRAL
C     (DERIVATIVES) FOR THE TOTALLY SYMMETRIC CASE
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      PARAMETER (LUTEMP = 48, IREPO = 0)
      LOGICAL LDIAG, FULMAT, LAEQB
      DIMENSION AO(*), SO(*)
      COMMON /CSYM1/ BUF(600), IBUF(600), LENGTH, INDMAX
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      COMMON /SYMIND/ INDFA(8,10), INDFB(8,10), ISOFRA(8), ISOFRB(8)
      COMMON /GENCON/ NRCA,NRCB, CONTA(MXCONT), CONTB(MXCONT)
      IBTAND(I,J) = AND(I,J)
      IBTXOR(I,J) = XOR(I,J)
      IF (IPRINT .GT. 10) THEN
       CALL HEADER('Subroutine SYM1S',-1)
       WRITE (LUPRI, '(A, I5)') ' KB ', KB
       WRITE (LUPRI, '(A,2I5)') ' NHKTA/B ', NHKTA, NHKTB
      END IF
C
C  LOOP OVER ALL IRREPS IN MOLECULE
C
      DO 100 IREP = 0, MAXLOP
       IF (FULMAT) THEN
        INDOFF = NPARSU(IREP + 1)
       ELSE
        IMAT = IMAT0 + IREP
       END IF
C
C  LOOP OVER AOS WHICH ARE OF SYMMETRY IREP IN STABILIZER MULA
C
       ISKIP=NRCA*NRCB
       IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
       DO 200 NA = 1, KHKTA
        IF (IBTAND(MULA,IBTXOR(IREP,ITYPE(NHKTA,NA))).EQ.0) THEN
         NAT    = KHKTB*(NA - 1)
         DO 300 NB = 1,KHKTB
          IF (IBTAND(MULB,IBTXOR(IREP,ITYPE(NHKTB,NB))).EQ.0) THEN
C
C              Weight and parity factor
C
          FAC = HKAB*PT(IBTAND(KB,IBTXOR(IREP,ITYPE(NHKTB,NB))))
C
C              Locate SO integrals to which AO's contribute
C
          INT=0
          DO 150 IRCA=1,NRCA
           MAXB=NRCB
           IF(LDIAG) MAXB=IRCA
           DO 150 IRCB=1,MAXB
            INT=INT+1
C
            IF(LAEQB.AND.IRCB.GT.IRCA) GO TO 150
            IF(IREP.EQ.0.AND.LDIAG.AND.(IRCA.EQ.IRCB).AND.(NA.LT.NB)
     *          ) GO TO 150
             IF(LAEQB.AND.(IRCA.EQ.IRCB).AND.NA.LT.NB) GO To 150
             INDA = INDFA(IREP + 1,NA)-NRCA+IRCA
             INDB = INDFB(IREP + 1,NB)-NRCB+IRCB
             INDM = MAX(INDA,INDB)
             FACT=FAC
             IF(INDA.LT.INDB) FACT=ASYM*FACT
             RINT = FACT*AO((INT+(NB+NAT-1)*ISKIP-1)*JSKIP+1)
             IF (FULMAT) THEN
              IND  = INDOFF + (INDM*(INDM - 3))/2 + INDA + INDB
              SO(IND) = SO(IND) + RINT
c              write(*,*)'ia,ib,value  ',inda,indb,so(ind),rint
             ELSE
              IF (ABS(RINT) .GT. THRESH) THEN
              IND  = (INDM*(INDM - 3))/2 + INDA + INDB
              INDMAX = MAX(IND,INDMAX)
              LABEL = IND*2**16 + IMAT
              IF (IPRINT .GT. 20) THEN
               WRITE (LUPRI,'(A,F12.6,2I3,2I2,I4,I2,I5)')
     *                'SYM1S - NA/B,IREPA/B,IND,IREPO,IMAT',
     *                 RINT, NA, NB, IREP, IREP, IND, IREPO,
     *                 IMAT
              END IF
              IF (LENGTH .EQ. 600) THEN
               WRITE (LUTEMP) BUF, IBUF, LENGTH
               IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     *            ' Buffer of length',LENGTH,
     *            ' has been written in SYM1S.'
               LENGTH = 0
              ENDIF
              LENGTH = LENGTH + 1
              BUF (LENGTH) = RINT
              IBUF(LENGTH) = LABEL
             ENDIF
            ENDIF
150        CONTINUE
          END IF
300      CONTINUE
        END IF
200    CONTINUE
100   CONTINUE
      RETURN
      END
