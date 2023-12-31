      SUBROUTINE G3AB1(GAMMA,BUF,GOFF,LENGAM,FIJKL)
C
C THIS ROUTINE PROCESSES THE G(Ij,Kl) GAMMA LISTS FOR SPIN CASE
C  AB AND FILLS THE SYMMETRY CASE ABAB MULLIKEN ORDERED GAMMA LIST.
C  THIS IS THE INCORE VERSION!
C
C NOTE THAT ABAB IN MULLIKEN CORRESPONDS TO AABB IN DIRAC.
C
C     G(IJ,kl) = G(IJ,kl) + G(JI,kl) + G(IJ,lk) + G(JI,lk)
C
C   THE GAMMAS ON DISK DO NOT HAVE THE COMPLETE SYMMETRY, BUT
C    ONLY SATISFY G(IJ,kl) = G(JI,lk).  CONSEQUENTLY, TWO
C    CONTRIBUTIONS TO EACH G(IJ,kl) ARE NEEDED.
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
C COMPUTE OFFSETS TO FIRST ELEMENT OF ALL G(XX,YY) IRREP POSITIONS
C
      LISTIJKL=113
      IOFFR=1
      IOFFL=1
      DO 1000 IRREPR=1,NIRREP
       IOFFSR(IRREPR)=IOFFR
       IOFFSL(IRREPR)=IOFFR
       IOFFR=IOFFR+POP(IRREPR,1)*POP(IRREPR,2)
1000  CONTINUE
      DO 100 IRREPDO=2,NIRREP
       DO 110 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPDO)
        IF(IRREPL.GT.IRREPR)THEN
         ILOGREC=IOFFSR(IRREPR)
         IFIRST =IOFFSL(IRREPL)
         NUML=POP(IRREPR,2)
         NUMK=POP(IRREPR,1)
         NUMJ=POP(IRREPL,2)
         NUMI=POP(IRREPL,1)
         MAXL=MOPOP(IRREPL)*MOPOP(IRREPR)
         MAXA=MOPOP(IRREPL)
         DO 120 INDXL=1,NUML
          DO 130 INDXK=1,NUMK
           INDL=INDXL
           INDK=INDXK
           CALL GETLST(BUF,ILOGREC,1,1,1,LISTIJKL)
           ILOGREC=ILOGREC+1
           IPOS=IFIRST
           DO 140 INDXJ=1,NUMJ
            DO 150 INDXI=1,NUMI
             X=BUF(IPOS)*FIJKL
             IPOS=IPOS+1
             INDI=INDXI
             INDJ=INDXJ
             GINDXD=INDXF(INDI,INDK,MAXA)
             GINDXR=INDXF(INDJ,INDL,MAXA)
             GINDX0=INDXF(GINDXD,GINDXR,MAXL)
             GINDX =GINDX0+GOFF(IRREPL,IRREPR)
             GAMMA(GINDX)=GAMMA(GINDX)+X
150         CONTINUE
140        CONTINUE
130       CONTINUE
120      CONTINUE
        ENDIF
110    CONTINUE
100   CONTINUE
C
C NOW WE MUST LOOP OVER ALL OTHER IRREPS AND PUT THE OTHER PIECE IN
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
       IOFFSL(IRREPR)=IOFFR
       IOFFR=IOFFR+POP(IRREPL,1)*POP(IRREPR,2)
1200  CONTINUE
      DO 170 IRREPB=1,NIRREP
       IRREPA=DIRPRD(IRREPB,IRREPDO)
       IF(IRREPA.LT.IRREPB)THEN
        ILOGREC=IOFFSR(IRREPB)
        IFIRST=IOFFSL(IRREPA)
        NUML=POP(IRREPB,2)
        NUMK=POP(IRREPA,1)
        NUMJ=POP(IRREPA,2)
        NUMI=POP(IRREPB,1)
        MAXL=MOPOP(IRREPA)*MOPOP(IRREPB)
        MAXA=MOPOP(IRREPA)
        MAXB=MOPOP(IRREPB)
        DO 180 INDXL=1,NUML
         DO 185 INDXK=1,NUMK
          INDL=INDXL
          INDK=INDXK
          CALL GETLST(BUF,ILOGREC,1,1,IRREPDO,LISTIJKL)
          ILOGREC=ILOGREC+1
          IPOS=IFIRST
           DO 190 INDXJ=1,NUMJ
            DO 195 INDXI=1,NUMI
             INDJ=INDXJ
             INDI=INDXI
             X=BUF(IPOS)*FIJKL
             IPOS=IPOS+1
             GINDXL=INDXF(INDI,INDK,MAXB)
             GINDXR=INDXF(INDL,INDJ,MAXB)
             GINDX0=INDXF(GINDXL,GINDXR,MAXL)
             GINDX=GINDX0+GOFF(IRREPB,IRREPA)
             GAMMA(GINDX)=GAMMA(GINDX)+X
195         CONTINUE
190        CONTINUE
185      CONTINUE
180     CONTINUE
       ENDIF
170   CONTINUE
160   CONTINUE
      RETURN
      END
