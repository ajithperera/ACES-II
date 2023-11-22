      SUBROUTINE PBLOCK(PSO,PAO,ICOMPA,ICOMPB,ICOMPC,ICOMPD,
     1                  ISOR,ISOS,ISOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     *********************************************************
C     ***  Transform a block of the P-matrix from SO basis  ***
C     ***  to AO basis.  This is done for fixed component   ***
C     ***  indices ICOMPx and fixed symmetry operations.    ***
C     *********************************************************
C
#include "mxcent.par"
#include "baslims.par"
      PARAMETER (ZERO = 0. D00)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      LOGICAL         AEQB,   CEQD,   DIAGAB, DIAGCD, DIACAC,
     *                ONECEN, PQSYM,  DTEST,
     *                TPRIAB, TPRICD, TCONAB, TCONCD
      COMMON /INTINF/ THRESH,
     *                NHKTA,  NHKTB,  NHKTC,  NHKTD,
     *                MAXAB,  MAXCD,  JMAX0,
     *                KHKTA,  KHKTB,  KHKTC,  KHKTD,
     *                KHKTAB, KHKTCD, KHABCD,
     *                MHKTA,  MHKTB,  MHKTC,  MHKTD,
     *                MULA,   MULB,   MULC,   MULD,
     *                NORBA,  NORBB,  NORBC,  NORBD, NORBAB, NORBCD,
     *                NUCA,   NUCB,   NUCC,   NUCD,  NUCAB,  NUCCD,
     *                NSETA,  NSETB,  NSETC,  NSETD,
     *                ISTEPA, ISTEPB, ISTEPC, ISTEPD,
     *                NSTRA,  NSTRB,  NSTRC,  NSTRD,
     *                AEQB,   CEQD,
     *                DIAGAB, IAB0X,  IAB0Y,  IAB0Z,
     *                DIAGCD, ICD0X,  ICD0Y,  ICD0Z,
     *                DIACAC, ONECEN, PQSYM,  IPQ0X, IPQ0Y, IPQ0Z,
     *                TPRIAB, TPRICD, TCONAB, TCONCD,
     *                MAXDER, DTEST
      DIMENSION PSO(1), PAO(1)
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
      IAR(I,J,K,L) = KHKTA*(KHKTB*(KHKTC*(L-1)+K-1)+J-1)+I
      IRR(I,J,K,L) = MULTA*(MULTB*(MULTC*(L)+K)+J)+I+1
      ISOTS = IBTXOR(ISOT,ISOS)
      MULTA = MULT(MULA)
      MULTB = MULT(MULB)
      MULTC = MULT(MULC)
      MULTD = MULT(MULD)
      NCABCD = NORBA*NORBB*NORBC*NORBD
c      MULTOT = IRR(MULTA-1,MULTB-1,MULTC-1,MULTD-1)
      MULTOT = IRR(MULTA-1,MULTB-1,MULTC-1,0)
      ITYNA = ITYPE(NHKTA,ICOMPA)
      ITYNB = ITYPE(NHKTB,ICOMPB)
      ITYNC = ITYPE(NHKTC,ICOMPC)
      ITYND = ITYPE(NHKTD,ICOMPD)
      DO 100 I = 1,NCABCD
         PAO(I) = ZERO
100   CONTINUE
   
CSSS      write(6, "(1x,2F20.10,1x,i4)") (PSO(i), i=1,NCABCD)
      

C
C     Loop over irreps
C
      IRCNTA = -1
      DO 10 IREPA = 0, MAXLOP
      IF(IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
         IRCNTA = IRCNTA + 1
         IRCNTB = -1
         DO 20 IREPB = 0, MAXLOP
         IF(IBTAND(MULB,IBTXOR(IREPB,ITYNB)) .EQ. 0) THEN
            IRCNTB = IRCNTB + 1
            IRCNTC = -1
            DO 30 IREPC = 0, MAXLOP
            IF(IBTAND(MULC,IBTXOR(IREPC,ITYNC)) .EQ. 0) THEN
               IRCNTC = IRCNTC + 1
               IRCNTD = -1
               DO 35 IREPD = 0,MAXLOP
               IF (IBTAND(MULD,IBTXOR(IREPD,ITYND)) .EQ. 0) THEN
                  IRCNTD = IRCNTD + 1
                  IF (IBTXOR(IREPD,IBTXOR(IREPA,IBTXOR(IREPB,IREPC)))
     1                   .NE. 0) GOTO 35
C
C                 Determine sign factor.  It is not necessary to put
C                 in factors from the angular part (factors like
C                 IBTXOR(IREPx,ITYNx) ) as these are folded into
C                 the integral derivatives already.
C
                  SIGN = PT(IBTAND(ISOR,IREPB))
     1                  *PT(IBTAND(ISOT,IREPC))
     2                  *PT(IBTAND(ISOTS,IREPD))
                  IOFF = NCABCD*(MULTOT
     1             *(IAR(ICOMPA,ICOMPB,ICOMPC,ICOMPD) - 1)
c     2             + IRR(IRCNTA,IRCNTB,IRCNTC,IRCNTD) - 1)
     2             + IRR(IRCNTA,IRCNTB,IRCNTC,0) - 1)
                  DO 40 I = 1,NCABCD
                      PAO(I) = PAO(I) + SIGN*PSO(IOFF+I)
CSSS                      Write(6,"(1x,2F20.10,1x,i4)"), PAO(I), PSO(IOFF+I)
CSSS     *                    ,ioff 
40                CONTINUE
               ENDIF
35          CONTINUE
            ENDIF
30          CONTINUE
         ENDIF
20       CONTINUE
      ENDIF
10    CONTINUE
CSSS      write(6, "(1x,2F20.10,1x,i4)") (PAO(i), i=1,NCABCD)
CSSS      Write(6,*)
      RETURN
      END
