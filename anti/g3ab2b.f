      SUBROUTINE G3AB2B(GAMMA,BUF,GOFF,LENGAM,FIJKA)
C
C THIS ROUTINE PROCESSES THE G(Ij,Ka) GAMMA LISTS FOR SPIN CASE
C  AB AND FILLS THE SYMMETRY CASE ABAB MULLIKEN ORDERED GAMMA LIST
C  FOR SPIN CASE AABB.  THIS CODE IS USED FOR UHF CALCULATIONS ONLY.
C
C NOTE THAT ABAB IN MULLIKEN CORRESPONDS TO AABB IN DIRAC.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION GAMMA(LENGAM),BUF(1),X,X1,X2 
      DOUBLE PRECISION FABIJ,FAIBJ,FIJKL,FABCD,FABCI,FIJKA
      INTEGER GOFF(8,8),IOFFSL(8),IOFFSR(8)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MOPOPS/ MOPOP(8)
C
      INDXF(I,J,N)=I+(J-1)*N
      INDXT(I,J)  =I+(J*(J-1))/2
      NNP1O2(I)   =(I*(I+1))/2
C
C FIRST DO LIST 110
C
C       21 21        22 11
C     G(IK,ja) <-- g<Ij,Ka>
C
C       21 21        11 22
C     G(KI,aj) <-- g<Ij,Ka>
C
C COMPUTE OFFSETS TO FIRST ELEMENT OF ALL G(XX,YY) IRREP POSITIONS
C  AND PUT IN FIRST CONTRIBUTION (THE PIECE FROM DPD IRREP = 1).
C
      LISTIJKA=110
      IOFFR=1
      IOFFL=1
      DO 1000 IRREPR=1,NIRREP
       IOFFSR(IRREPR)=IOFFR
       IOFFSL(IRREPR)=IOFFL
       IOFFR=IOFFR+POP(IRREPR,1)*VRT(IRREPR,2)
       IOFFL=IOFFL+POP(IRREPR,1)*POP(IRREPR,2)
1000  CONTINUE
      DO 100 IRREPDO=2,NIRREP
       DO 110 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPDO)
        ILOGREC=IOFFSR(IRREPR)
        IFIRST =IOFFSL(IRREPL)
        NUMA=VRT(IRREPR,2)
        NUMK=POP(IRREPR,1)
        NUMJ=POP(IRREPL,2)
        NUMI=POP(IRREPL,1)
        MAXL=MOPOP(IRREPL)*MOPOP(IRREPR)
        MAXA=MOPOP(IRREPL)
        MAXB=MOPOP(IRREPR)
        DO 120 INDXA=1,NUMA
         DO 130 INDXK=1,NUMK
          INDA=INDXA+POP(IRREPR,2)
          INDK=INDXK
          CALL GETLST(BUF,ILOGREC,1,1,1,LISTIJKA)
          ILOGREC=ILOGREC+1
          IPOS=IFIRST
          IF(IRREPL.GT.IRREPR)THEN
C
C       21 21        22 11
C     G(IK,ja) <-- g<Ij,Ka>
C
           DO 140 INDXJ=1,NUMJ
            DO 150 INDXI=1,NUMI
             X=BUF(IPOS)*FIJKA
             IPOS=IPOS+1
             INDI=INDXI
             INDJ=INDXJ
             GINDXL=INDXF(INDI,INDK,MAXA)
             GINDXR=INDXF(INDJ,INDA,MAXA)
             GINDX0=INDXF(GINDXR,GINDXL,MAXL)
             GINDX =GINDX0+GOFF(IRREPL,IRREPR)
             GAMMA(GINDX)=GAMMA(GINDX)+X
150         CONTINUE
140        CONTINUE
          ELSE
C
C       21 21        11 22
C     G(KI,aj) <-- g<Ij,Ka>
C
           IPOS=IFIRST
           DO 141 INDXJ=1,NUMJ
            DO 151 INDXI=1,NUMI
             X=BUF(IPOS)*FIJKA
             IPOS=IPOS+1
             INDI=INDXI
             INDJ=INDXJ
             GINDXL=INDXF(INDK,INDI,MAXB)
             GINDXR=INDXF(INDA,INDJ,MAXB)
             GINDX0=INDXF(GINDXR,GINDXL,MAXL)
             GINDX =GINDX0+GOFF(IRREPR,IRREPL)
             GAMMA(GINDX)=GAMMA(GINDX)+X
151         CONTINUE
141        CONTINUE
          ENDIF
130      CONTINUE
120     CONTINUE
110    CONTINUE
100   CONTINUE
C
C NOW WE MUST LOOP OVER ALL OTHER IRREPS AND PUT THE OTHER PIECES IN
C
C       21 21        21 12
C     G(IK,aj) <-- g<Ij,Ka>
C
C       21 21        12 21
C     G(KI,ja) <-- g<Ij,Ka>
C
C
      DO 160 IRREPDO=2,NIRREP
C
C COMPUTE OFFSETS TO ALL BA AB POSITIONS
C
      IOFFR=1
      IOFFL=1
      DO 1200 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPDO)
       IOFFSR(IRREPR)=IOFFR
       IOFFSL(IRREPR)=IOFFL
       IOFFR=IOFFR+POP(IRREPL,1)*VRT(IRREPR,2)
       IOFFL=IOFFL+POP(IRREPL,1)*POP(IRREPR,2)
1200  CONTINUE
      DO 170 IRREPB=1,NIRREP
       IRREPA=DIRPRD(IRREPB,IRREPDO)
       ILOGREC=IOFFSR(IRREPB)
       IFIRST=IOFFSL(IRREPA)
       NUMA=VRT(IRREPB,2)
       NUMK=POP(IRREPA,1)
       NUMJ=POP(IRREPA,2)
       NUMI=POP(IRREPB,1)
       MAXL=MOPOP(IRREPA)*MOPOP(IRREPB)
       MAXA=MOPOP(IRREPA)
       MAXB=MOPOP(IRREPB)
       DO 180 INDXA=1,NUMA
        DO 185 INDXK=1,NUMK
         INDA=INDXA+POP(IRREPB,2)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPDO,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IFIRST
         IF(IRREPA.LT.IRREPB)THEN
C
C       21 21        21 12
C     G(IK,aj) <-- g<Ij,Ka>
C
          DO 190 INDXJ=1,NUMJ
           DO 195 INDXI=1,NUMI
            INDJ=INDXJ
            INDI=INDXI
            X=BUF(IPOS)*FIJKA
            IPOS=IPOS+1
            GINDXL=INDXF(INDI,INDK,MAXB)
            GINDXR=INDXF(INDA,INDJ,MAXB)
            GINDX0=INDXF(GINDXR,GINDXL,MAXL)
            GINDX=GINDX0+GOFF(IRREPB,IRREPA)
            GAMMA(GINDX)=GAMMA(GINDX)+X
195        CONTINUE
190       CONTINUE
         ELSE
C
C       21 21        12 21
C     G(KI,ja) <-- g<Ij,Ka>
C
          DO 191 INDXJ=1,NUMJ
           DO 196 INDXI=1,NUMI
            INDJ=INDXJ
            INDI=INDXI
            X=BUF(IPOS)*FIJKA
            IPOS=IPOS+1
            GINDXL=INDXF(INDK,INDI,MAXA)
            GINDXR=INDXF(INDJ,INDA,MAXA)
            GINDX0=INDXF(GINDXR,GINDXL,MAXL)
            GINDX=GINDX0+GOFF(IRREPA,IRREPB)
            GAMMA(GINDX)=GAMMA(GINDX)+X
196        CONTINUE
191       CONTINUE
         ENDIF
185     CONTINUE
180    CONTINUE
170   CONTINUE
160   CONTINUE
C
C NOW DO LIST 109
C
C       21 21        22 11
C     G(IA,jk) <-- g<Ij,Ak>
C
C       21 21        11 22
C     G(AI,kj) <-- g<Ij,Ak>
C
C COMPUTE OFFSETS TO FIRST ELEMENT OF ALL G(XX,YY) IRREP POSITIONS
C  AND PUT IN FIRST CONTRIBUTION (THE PIECE FROM DPD IRREP = 1).
C
      LISTIJKA=109
      IOFFR=1
      IOFFL=1
      DO 2000 IRREPR=1,NIRREP
       IOFFSR(IRREPR)=IOFFR
       IOFFSL(IRREPR)=IOFFL
       IOFFR=IOFFR+VRT(IRREPR,1)*POP(IRREPR,2)
       IOFFL=IOFFL+POP(IRREPR,1)*POP(IRREPR,2)
2000  CONTINUE
      DO 200 IRREPDO=2,NIRREP
       DO 210 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPDO)
        ILOGREC=IOFFSR(IRREPR)
        IFIRST =IOFFSL(IRREPL)
        NUMK=POP(IRREPR,2)
        NUMA=VRT(IRREPR,1)
        NUMJ=POP(IRREPL,2)
        NUMI=POP(IRREPL,1)
        MAXL=MOPOP(IRREPL)*MOPOP(IRREPR)
        MAXA=MOPOP(IRREPL)
        MAXB=MOPOP(IRREPR)
        DO 220 INDXK=1,NUMK
         DO 230 INDXA=1,NUMA
          INDA=INDXA+POP(IRREPR,1)
          INDK=INDXK
          CALL GETLST(BUF,ILOGREC,1,1,1,LISTIJKA)
          ILOGREC=ILOGREC+1
          IPOS=IFIRST
          IF(IRREPL.GT.IRREPR)THEN
C
C       21 21        22 11
C     G(IA,jk) <-- g<Ij,Ak>
C
           DO 240 INDXJ=1,NUMJ
            DO 250 INDXI=1,NUMI
             X=BUF(IPOS)*FIJKA
             IPOS=IPOS+1
             INDI=INDXI
             INDJ=INDXJ
             GINDXL=INDXF(INDI,INDA,MAXA)
             GINDXR=INDXF(INDJ,INDK,MAXA)
             GINDX0=INDXF(GINDXR,GINDXL,MAXL)
             GINDX =GINDX0+GOFF(IRREPL,IRREPR)
             GAMMA(GINDX)=GAMMA(GINDX)+X
250         CONTINUE
240        CONTINUE
          ELSE
C
C       21 21        11 22
C     G(AI,kj) <-- g<Ij,Ak>
C
           DO 241 INDXJ=1,NUMJ
            DO 251 INDXI=1,NUMI
             X=BUF(IPOS)*FIJKA
             IPOS=IPOS+1
             INDI=INDXI
             INDJ=INDXJ
             GINDXL=INDXF(INDA,INDI,MAXB)
             GINDXR=INDXF(INDK,INDJ,MAXB)
             GINDX0=INDXF(GINDXR,GINDXL,MAXL)
             GINDX =GINDX0+GOFF(IRREPR,IRREPL)
             GAMMA(GINDX)=GAMMA(GINDX)+X
251         CONTINUE
241        CONTINUE
          ENDIF
230      CONTINUE
220     CONTINUE
210    CONTINUE
200   CONTINUE
C
C NOW WE MUST LOOP OVER ALL OTHER IRREPS AND PUT THE OTHER PIECES IN
C
C       21 21        21 12
C     G(IA,kj) <-- g<Ij,Ak>
C
C       21 21        12 21
C     G(AI,jk) <-- g<Ij,Ak>
C
C
      DO 260 IRREPDO=2,NIRREP
C
C COMPUTE OFFSETS TO ALL BA AB POSITIONS
C
      IOFFR=1
      IOFFL=1
      DO 2200 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPDO)
       IOFFSR(IRREPR)=IOFFR
       IOFFSL(IRREPR)=IOFFL
       IOFFR=IOFFR+VRT(IRREPL,1)*POP(IRREPR,2)
       IOFFL=IOFFL+POP(IRREPL,1)*POP(IRREPR,2)
2200  CONTINUE
      DO 270 IRREPB=1,NIRREP
       IRREPA=DIRPRD(IRREPB,IRREPDO)
       ILOGREC=IOFFSR(IRREPB)
       IFIRST=IOFFSL(IRREPA)
       NUMK=POP(IRREPB,2)
       NUMA=VRT(IRREPA,1)
       NUMJ=POP(IRREPA,2)
       NUMI=POP(IRREPB,1)
       MAXL=MOPOP(IRREPA)*MOPOP(IRREPB)
       MAXA=MOPOP(IRREPA)
       MAXB=MOPOP(IRREPB)
       DO 280 INDXK=1,NUMK
        DO 285 INDXA=1,NUMA
         INDA=INDXA+POP(IRREPA,1)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPDO,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IFIRST
         IF(IRREPA.LT.IRREPB)THEN
C
C       21 21        21 12
C     G(IA,kj) <-- g<Ij,Ak>
C
          DO 290 INDXJ=1,NUMJ
           DO 295 INDXI=1,NUMI
            INDJ=INDXJ
            INDI=INDXI
            X=BUF(IPOS)*FIJKA
            IPOS=IPOS+1
            GINDXL=INDXF(INDI,INDA,MAXB)
            GINDXR=INDXF(INDK,INDJ,MAXB)
            GINDX0=INDXF(GINDXR,GINDXL,MAXL)
            GINDX=GINDX0+GOFF(IRREPB,IRREPA)
            GAMMA(GINDX)=GAMMA(GINDX)+X
295        CONTINUE
290       CONTINUE
         ELSE
C
C       21 21        12 21
C     G(AI,jk) <-- g<Ij,Ak>
C
          DO 291 INDXJ=1,NUMJ
           DO 296 INDXI=1,NUMI
            INDJ=INDXJ
            INDI=INDXI
            X=BUF(IPOS)*FIJKA
            IPOS=IPOS+1
            GINDXL=INDXF(INDA,INDI,MAXA)
            GINDXR=INDXF(INDJ,INDK,MAXA)
            GINDX0=INDXF(GINDXR,GINDXL,MAXL)
            GINDX=GINDX0+GOFF(IRREPA,IRREPB)
            GAMMA(GINDX)=GAMMA(GINDX)+X
296        CONTINUE
291       CONTINUE
         ENDIF
285     CONTINUE
280    CONTINUE
270   CONTINUE
260   CONTINUE
      RETURN
      END
