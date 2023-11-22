
        


      SUBROUTINE SYMOUT(SO,NOABCD,SHAEQB,SHCEQD,SHABAB,FIRST,LAST,
     *                  IPRINT)
C
C     Write out blocks of symmetry integrals, eliminating duplicates
C
C                                          880412   PRT
C
C     Some low-brain work has been done by TUH
C
C     Rewritten to allow for triangular looping and to eliminate
C     all duplicates by testing  880601 TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "mxcent.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      DIMENSION SO(1), BUF(600), IBUF(600)
      LOGICAL SHAEQB, SHCEQD, SHABAB,
     *        DCMPAB, DCMPCD, DCMPAC,
     *        DRALTB, DRCLTD, FIRST, LAST,
     *        DRABAB, DCABAB, IAEQIC, IALTIC
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
      COMMON /HRUNIT/ LUONEL, LUINTA, LUSUPM, LUDASP, LUSOL
      SAVE BUF, IBUFI, ICOUNT, NBUF
C
      CHARACTER*8 BASLBL(4)
      DATA BASLBL /'********','        ','        ','BASTWOEL'/
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (IPRINT .GT. 6) THEN
         CALL HEADER('Subroutine SYMOUT',-1)
      END IF
C
C     *************************************
C     ***** Last call to empty buffer *****
C     *************************************
C
      IF (LAST) THEN
         IF (ICOUNT .GT. 0) THEN
            WRITE (LUINTA) BUF, IBUF, ICOUNT
         END IF
         WRITE (LUINTA) BUF, IBUF, -1
         REWIND LUINTA
         CLOSE (UNIT=LUINTA)
         NTOTAL = 600*NBUF + ICOUNT
         WRITE (LUPRI,'(/1X,I10,A//)') NTOTAL,
     *      ' undifferentiated two-electron integrals'//
     *      ' have been written on disk.'
         RETURN
      END IF
C
C     *******************************************************
C     ***** Initialization when subroutine first called *****
C     *******************************************************
C
      IF (FIRST) THEN
         REWIND LUINTA
         WRITE (LUINTA) BASLBL
         ICOUNT = 0
         NBUF = 0
      END IF
C
C     **********************************************************
C     ***** Loop over shell components for this quadruplet *****
C     **********************************************************
C
      IF (IPRINT .GT. 6) THEN
         WRITE(LUPRI,'(/4X,A,6X,A,5X,A,6X,A,3X,A/)')
     *     'Components', 'Irreps', 'Sequence', 'Indices', 'Integral'
      END IF
      ISOFF = 0
      DO 100 NA = 1,KHKTA
        NSTRNA = NSTRA + NA
        ITYNA = ITYPE(NHKTA,NA)
        KHKTBB = KHKTB
        IF (DIAGAB) KHKTBB = NA
        DO 110 NB = 1,KHKTBB
          NSTRNB = NSTRB + NB
          ITYNB = ITYPE(NHKTB,NB)
          DCMPAB = SHAEQB .AND. NA .EQ. NB
          KHKTCC = KHKTC
          IF (SHABAB) KHKTCC = NA
          DO 120 NC = 1,KHKTCC
            NSTRNC = NSTRC + NC
            ITYNC = ITYPE(NHKTC,NC)
            DCMPAC = SHABAB .AND. NA .EQ. NC
            KHKTDD = KHKTD
            IF (DIAGCD) KHKTDD = NC
            IF (DCMPAC) KHKTDD = NB
            DO 130 ND = 1,KHKTDD
              NSTRND = NSTRD + ND
              ITYND = ITYPE(NHKTD,ND)
              DCMPCD = SHCEQD .AND. NC .EQ. ND
              DCABAB = DCMPAC .AND. NB .EQ. ND
C
C             Loop over all irreps in molecule
C
              DO 200 IREPA = 0, MAXLOP
              IF (IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
                DO 210 IREPB = 0, MAXLOP
                IF (IBTAND(MULB,IBTXOR(IREPB,ITYNB)) .EQ. 0) THEN
                  DRALTB = DCMPAB .AND. IREPA .LT. IREPB
                  DO 220 IREPC = 0, MAXLOP
                  IF (IBTAND(MULC,IBTXOR(IREPC,ITYNC)) .EQ. 0) THEN
                    IREPD = IBTXOR(IBTXOR(IREPA,IREPB),IREPC)
                  IF (IBTAND(MULD,IBTXOR(IREPD,ITYND)) .EQ. 0) THEN
                  IF (DCABAB .AND. (IREPA .LT. IREPC .OR.
     *               (IREPA.EQ.IREPC .AND. IREPB.LT.IREPD))) GOTO 230
                    DRCLTD = DCMPCD .AND. IREPC .LT. IREPD
                    DRABAB = DCABAB .AND. IREPA .EQ. IREPC
     *                              .AND. IREPB .EQ. IREPD
C
C                   *****************************************
C                   ***** Loop over contracted orbitals *****
C                   *****************************************
C
                    INT = 0
                    DO 300 IA = 1,NORBA
                      INDA = IPTSYM(NSTRNA + KHKTA*(IA-1),IREPA)
                      NORBBB = NORBB
                      IF (TCONAB) NORBBB = IA
                    DO 310 IB = 1,NORBBB
                      IF (DCMPAB) THEN
                         IF ((IB.GT.IA) .OR. (DRALTB.AND.IB.EQ.IA)) THEN
                            INT = INT + NORBCD
                            GO TO 310
                         END IF
                      END IF
                      INDB = IPTSYM(NSTRNB + KHKTB*(IB-1),IREPB)
                      INDAB = MAX(INDA,INDB)*255 + INDA + INDB
                    DO 320 IC = 1,NORBC
                      INDC = IPTSYM(NSTRNC + KHKTC*(IC-1),IREPC)
                      IAEQIC = DRABAB .AND. IA .EQ. IC
                      IALTIC = DRABAB .AND. IA .LT. IC
                      NORBDD = NORBD
                      IF (TCONCD) NORBDD = IC
                    DO 330 ID = 1,NORBDD
                      INT = INT + 1
                      IF (DCMPCD ) THEN
                         IF (ID.GT.IC) GO TO 330
                         IF (DRCLTD .AND. ID.EQ.IC) GO TO 330
                      END IF
                      IF ((IALTIC .OR. (IAEQIC.AND.IB.LT.ID))) GOTO 330
                      SOINT = SO(ISOFF+INT)
                      IF (ABS(SOINT) .GT. THRESH) THEN
                         INDD = IPTSYM(NSTRND + KHKTD*(ID-1),IREPD)
                         IF (IPRINT .GT. 6) THEN
                           WRITE(LUPRI,
     *                     '(4I3,2X,4I3,2X,4I3,2X,4I3,2X,1P,D16.8)')
     *                        NA,    NB,    NC,    ND,
     *                        IREPA, IREPB, IREPC, IREPD,
     *                        IA,    IB,    IC,    ID,
     *                        INDA,  INDB,  INDC,  INDD, SOINT
                         END IF
                         INDCD  = MAX(INDC,INDD)*255 + INDC + INDD
                         IABCD  = MAX(INDAB,INDCD)*65535 + INDAB + INDCD
                         ICOUNT = ICOUNT + 1
                         BUF (ICOUNT) = SOINT
                         IBUF(ICOUNT) = IABCD
                         IF (ICOUNT .EQ. 600) THEN
                            WRITE (LUINTA) BUF,IBUF,ICOUNT
                            NBUF = NBUF + 1
                            ICOUNT = 0
                            IF (IPRINT .GE. 6) THEN
                                WRITE (LUPRI,'(2X,A,I5,A/)')
     *                                'Integral buffer #',NBUF,
     *                                'has been written.'
                            END IF
                         END IF
                      END IF
330                 CONTINUE
320                 CONTINUE
310                 CONTINUE
300                 CONTINUE
230                 CONTINUE
                    ISOFF = ISOFF + NOABCD
                  END IF
                  END IF
220               CONTINUE
                END IF
210             CONTINUE
              END IF
200           CONTINUE
130         CONTINUE
120       CONTINUE
110     CONTINUE
100   CONTINUE
      RETURN
      END
