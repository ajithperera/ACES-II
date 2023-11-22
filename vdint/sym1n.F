      SUBROUTINE SYM1N(AO,SO,IREPO,KB,MULA,MULB,NHKTA,NHKTB,
     &                 KHKTA,KHKTB,HKAB,LDIAG,LAEQB,FULMAT,
     &                 THRESH,IMAT0,JSKIP,IPRINT,ASYM)
C
C  Take block of distinct AO integral (derivatives) and
C  generate symmetrized contributions to SO integral
C  (derivatives) over non-symmetric operators
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL LDIAG, FULMAT, LAEQB
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      PARAMETER (LUTEMP = 48)
C
      DIMENSION AO(*), SO(*)
      COMMON /CSYM1/ BUF(600), IBUF(600), LENGTH, INDMAX
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
      COMMON /SYMIND/ INDFA(8,10), INDFB(8,10), ISOFRA(8), ISOFRB(8)
      COMMON /GENCON/ NRCA,NRCB,CONTA(MXCONT), CONTB(MXCONT)
C
      IBTAND(I,J) = AND(I,J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GT. 10) THEN
       CALL HEADER('Subroutine SYM1N',-1)
       WRITE (LUPRI, '(A,2I5)') ' IREPO, KB ', IREPO, KB
       WRITE (LUPRI, '(A,2I5)') ' NHKTA/B ', NHKTA, NHKTB
      END IF
C
C  Loop over irreps for first basis function - those for second
C  are obtained from operator symmetry IREPO
C
      DO 100 IREPA = 0, MAXLOP
       IREPB = IBTXOR(IREPO,IREPA)
       IF (FULMAT) INDOFF = NPARNU(IREPO+1,MAX(IREPA,IREPB)+1)
C
C  Loop over AOs which are of symmetry IREPA in stabilizer MULA
C
       ISKIP=NRCA*NRCB
       IF(LDIAG) ISKIP=NRCA*(NRCA+1)/2
C
       DO 200 NA = 1, KHKTA
        IF (IBTAND(MULA,IBTXOR(IREPA,ITYPE(NHKTA,NA))).EQ.0) THEN
         NAT = KHKTB*(NA - 1)
         DO 300 NB = 1,KHKTB
          IF (IBTAND(MULB,IBTXOR(IREPB,ITYPE(NHKTB,NB))).EQ.0) THEN
C
C  Weight and parity factor
C
           FAC = HKAB*PT(IBTAND(KB,IBTXOR(IREPB,ITYPE(NHKTB,NB))))
C
C  Locate SO integrals to which AOs contribute
C
           INT=0
           DO 160 IRCA=1,NRCA
            MAXB=NRCB
            IF(LDIAG) MAXB=IRCA
            DO 150 IRCB=1,MAXB
             INT=INT+1
             IF(LDIAG.AND.(IRCA.EQ.IRCB).AND.NA.LT.NB) GO TO 150
             IF(LAEQB.AND.NA.EQ.NB.AND.IRCA.EQ.IRCB.AND.
     &          IREPA.LT.IREPB) GO TO  150
             IF(LAEQB.AND.IRCA.LT.IRCB) GO TO 150
             IF(LAEQB.AND.(IRCA.EQ.IRCB).AND.NA.LT.NB) GO TO 150
             INDA = INDFA(IREPA + 1,NA)-NRCA+IRCA
             INDB = INDFB(IREPB + 1,NB)-NRCB+IRCB
             RINT = FAC*AO((INT+(NAT+NB-1)*ISKIP-1)*JSKIP+1)
             IF (FULMAT) THEN
              IF (IREPA .GE. IREPB) THEN
               IND  = INDOFF + NAOS(IREPB+1)*(INDA-1) + INDB
               RINT=RINT*ASYM
              ELSE
               IND  = INDOFF + NAOS(IREPA+1)*(INDB-1) + INDA
              ENDIF
              SO(IND) = SO(IND) + RINT
             ELSE
              IF (ABS(RINT) .GT. THRESH) THEN
               IF (IREPA .GE. IREPB) THEN
                IND  = NAOS(IREPB+1)*(INDA-1) + INDB
               ELSE
                IND  = NAOS(IREPA+1)*(INDB-1) + INDA
               ENDIF
               INDMAX = MAX(IND,INDMAX)
               LABEL  = IND*2**16 + IMAT0 + MAX(IREPA,IREPB)
               IF (IPRINT .GT. 20) THEN
                WRITE (LUPRI,'(A,F12.6,2I3,2I2,I4,I2,I5)')
     &                 'SYM1N - NA/B,IREPA/B,IND,IREPO,IMAT',
     &                 RINT, NA, NB, IREPA, IREPB, IND, IREPO,
     &                 IMAT0 + MAX(IREPA,IREPB)
               END IF
               IF (LENGTH .EQ. 600) THEN
                WRITE (LUTEMP) BUF, IBUF, LENGTH
                IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     &                             ' Buffer of length',LENGTH,
     &                             ' has been written in SYM1N.'
                LENGTH = 0
               ENDIF
               LENGTH = LENGTH + 1
               BUF(LENGTH)  = RINT
               IBUF(LENGTH) = LABEL
              ENDIF
             ENDIF
150         CONTINUE
160        CONTINUE
          END IF
300      CONTINUE
        END IF
200    CONTINUE
100   CONTINUE
      RETURN
      END
