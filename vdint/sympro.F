      SUBROUTINE SYMPRO(CTRAN,ITRAN,JPRX,JTRAN,KWO,
     *                  MAMN,NB,NORB,NSYM,LAMN,IPRINT)
C
C     This subroutine determines the following parameters:
C
C     CTRAN, IPIND, IPTSYM, ITRAN, JPRX,
C     JTRAN, MAMN, NAOS, NB, NPAR, NPARNU, NPARSU, NORB, NSYM
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      LOGICAL SHARE
      LOGICAL DCORD, DCORGD, NOORBT, DOPERT
      CHARACTER NAMEX*6
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "mxcent.par"
#include "baslims.par"
C
      CHARACTER*4 KWO(20), JPRX(MXCORB), LAMN(MXSHEL), MAMN(MXCORB)
C
      DIMENSION JTRAN(MXCORB), CTRAN(MXCORB,8), ITRAN(MXCORB,8), NB(8),
     &          IRREP(MXCORB)
C
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
      COMMON /SHELLSi/ KMAX,
     &                NHKT(MXSHEL),   KHKT(MXSHEL), MHKT(MXSHEL),
     &                ISTBAO(MXSHEL), NUCO(MXSHEL), JSTRT(MXSHEL),
     &                NSTRT(MXSHEL),  MST(MXSHEL),  NCENT(MXSHEL),
     &                NRCO(MXSHEL), NUMCF(MXSHEL),
     &                NBCH(MXSHEL),   KSTRT(MXSHEL)
      COMMON /SHELLS/ CENT(MXSHEL,3), SHARE(MXSHEL)
      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(8,36), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXCENT,0:7), NCRREP(0:7),
     &                IPTCOR(MXCENT*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
      COMMON /PINCOM/ IPIND(MXCORB), IBLOCK(MXCORB), INDGEN(MXCORB)
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      NORB      = 0
      NPARSU(1) = 0
      NSYM      = 0
C
C Run over Irreducible Representations LAMBDA
C
      DO 100 LA = 0, MAXLOP
       KABAS    = 0
       NUMR     = 0
       NAOS(LA+1) = 0
C
C Run over Symmetry Independent Shells
C
       NORBD  = 0
       ISHELL = 0
       DO 110 IA = 1, KMAX
        MULA  = ISTBAO(IA)
        NHKTA = NHKT(IA)
        KHKTA = KHKT(IA)
        MHKTA = MHKT(IA)
        NRCA  = NRCO(IA)
        NAA   = 0
        NAAA  = 0
        NUMCFA = NUMCF(IA)
        IF (NUMCFA .EQ. 1) ISHELL = ISHELL + 1
C
C run over all types of primitives within this shell
C
         DO 120 NA = 1, NHKTA*(NHKTA + 1)/2
          NAA  = NAA  + 1
C
C run over all contractions for this shell
C
          DO 121 IRCA = 1,NRCA
           IF (IBTAND(IBTSHR(MHKTA,NA-1),1) .EQ. 1) THEN
            IVARB = IBTXOR(LA,ITYPE(NHKTA,NA))
            NAAA = NAAA + 1
            NUMR = NUMR + 1
            NORBD = NORBD + 1
            IF (LA.EQ.0) IRREP(NUMR) = 0
            JKB = 0
C
C loop over all irreducible representations 
C
            DO 130 KB = 0, MAXLOP
             IF (IBTAND(KB,MULA) .EQ. 0) THEN
              JKB = JKB + 1
              KABAS = KABAS + 1
              IF (IBTAND(MULA,IVARB) .EQ. 0) THEN
               IF (NORB.LE.MXCORB) THEN
                CTRAN(NORB+1,JKB) = PT(IBTAND(KB,IVARB))
                ITRAN(NORB+1,JKB) = KABAS
               END IF
              END IF
             END IF
130         CONTINUE
C
C If orbital contributes to this representation:
C
            IF (IBTAND(MULA,IVARB) .EQ. 0) THEN
             NORB = NORB + 1
             NAOS(LA+1) = NAOS(LA+1) + 1
             MXYZ = ((NHKTA+1)*NHKTA*(NHKTA-1))/6 + NA
             JPRX(NORB) = KWO(MXYZ)
             MAMN(NORB) = LAMN(IA)
             IRREPN = IRREP(NUMR)
             IRREP(NUMR) = IRREPN + 1
             IPIND(NORB) = IA*2**16 + NAAA*2**8 + IRREPN
             INDGEN(NORB) = IRREPN*2**24 + NAA*2**16
     &                    + IRCA*2**8 + ISHELL
             JTRAN(NORB) = JKB
             IPTSYM(NORBD,LA) = NORB
            ELSE
             IPTSYM(NORBD,LA) = 0
            END IF
           END IF
121       CONTINUE
120      CONTINUE
110     CONTINUE
        IF (LA .GT. 0) NPARSU(LA+1) = NPARSU(LA) + NPAR(LA)
        NBASI      = NAOS(LA+1)
        NPAR(LA+1) = NBASI*(NBASI + 1)/2
        IF (NBASI .GT. 0) THEN
         NSYM     = NSYM + 1
         NB(NSYM) = NBASI
        END IF
100    CONTINUE
C
C Compute symmetry offsets for non-symmetric operators
C (Note case IREPO = 0 is thereby excluded)
C
       DO 200 IREPO = 1,MAXLOP
        IBLK = 0
        DO 210 IREPA = 0,MAXLOP
         IREPB = IBTXOR(IREPO,IREPA)
         IF (IREPA .GT. IREPB) THEN
          NPARNU(IREPO+1,IREPA+1) = IBLK
          IBLK = IBLK + NAOS(IREPA+1)*NAOS(IREPB+1)
         ENDIF
210     CONTINUE
200    CONTINUE
      RETURN
      END
