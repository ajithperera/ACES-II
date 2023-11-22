      SUBROUTINE SYM2(SO,AO,ISYMR,ISYMT,ISYMTS,SHABAB,NOABCD,FACTOR,
     *                IREPE,ISOFF,IPRINT)
C
C     Take block of distinct AO two-electron integral derivatives and
C     generate symmetrized contributions to SO integral derivatives
C
C                                          880429   PRT & TUH
C
C     Modified for triangular looping in components tuh
C     Loop over operator irreps outside loop over components 880928 tuh
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      LOGICAL SHABAB
      DIMENSION AO(*), SO(*)
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
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     *                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     *                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     *                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     *                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     *                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     *                IPTXYZ(3,0:7)
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     Loop over components
C
      IAOFF = 0
      ISOFF = 0
      DO 100 ICOMPA = 1,KHKTA
        ITYNA = ITYPE(NHKTA,ICOMPA)
        KHKTBB = KHKTB
        IF (DIAGAB) KHKTBB = ICOMPA
      DO 200 ICOMPB = 1,KHKTBB
        ITYNB = ITYPE(NHKTB,ICOMPB)
      DO 300 ICOMPC = 1,KHKTC
        ITYNC = ITYPE(NHKTC,ICOMPC)
        KHKTDD = KHKTD
        IF (DIAGCD) KHKTDD = ICOMPC
      DO 400 ICOMPD = 1,KHKTDD
        ITYND = ITYPE(NHKTD,ICOMPD)
        IF (SHABAB .AND. (ICOMPA .LT. ICOMPC .OR.
     *       (ICOMPA .EQ. ICOMPC .AND. ICOMPB .LT. ICOMPD))) GO TO 240
C
C       Loop over irreps
C
        DO 500 IREPA = 0, MAXLOP
          IF (IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .NE. 0) GOTO 500
          IRPAE = IBTXOR(IREPA,IREPE)
          DO 600 IREPB = 0, MAXLOP
            IF (IBTAND(MULB,IBTXOR(IREPB,ITYNB)) .NE. 0) GOTO 600
            IRPABE = IBTXOR(IREPB,IRPAE)
            FACB   = FACTOR*PT(IBTAND(ISYMR,IBTXOR(IREPB,ITYNB)))
            DO 700 IREPC = 0, MAXLOP
              IF (IBTAND(MULC,IBTXOR(IREPC,ITYNC)) .NE. 0) GOTO 700
              IREPD = IBTXOR(IREPC,IRPABE)
              IF (IBTAND(MULD,IBTXOR(IREPD,ITYND)) .NE. 0) GOTO 700
              FAC=PT(IBTXOR(IBTAND(ISYMT, IBTXOR(IREPC,ITYNC)),
     *                      IBTAND(ISYMTS,IBTXOR(IREPD,ITYND))))*FACB
C             IF (IPRINT .GT. 25) THEN
C               WRITE (LUPRI,'(A,4I5,5x,4I5)') ' comps, IREPE/A/B/C ',
C    *                  ICOMPA, ICOMPB, ICOMPC, ICOMPD,
C    *                  IREPE, IREPA, IREPB, IREPC
C               WRITE(LUPRI,'(A,6F8.4)') 'FACTOR, FAC', FACTOR, FAC
C             END IF
C
C             Compute SO integral contributions from this block
C             of contracted AO integrals.
C
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
              DO 800 INT = 1,NOABCD
                SO(ISOFF+INT) = SO(ISOFF+INT) + FAC*AO(IAOFF+INT)
800           CONTINUE
              ISOFF = ISOFF + NOABCD
700         CONTINUE
600       CONTINUE
500     CONTINUE
240   CONTINUE
      IAOFF = IAOFF + NOABCD
400   CONTINUE
300   CONTINUE
200   CONTINUE
100   CONTINUE
      RETURN
      END
