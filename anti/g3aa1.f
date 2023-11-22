      SUBROUTINE G3AA1(GAMMA,BUF,GOFF,LENGAM,ISPIN,FIJKL)
C
C THIS ROUTINE PROCESSES THE G(I<J,K<L) GAMMA LISTS FOR SPIN CASES
C  AAAA AND BBBB AND FILLS THE SYMMETRY CASE ABAB MULLIKEN ORDERED 
C  GAMMA LIST.  THIS CODE IS USED FOR UHF CALCULATIONS ONLY.
C
C       21 21        22 11
C     G(IK,JL) <-- g<IJ,KL>
C  
C       21 21        12 21        21 21
C     G(KI,JL) <-- g<IJ,KL> = - g<JI,KL>
C   
C IN ADDITION, THE FOLLOWING SYMMETRIZATION IS REQUIRED:
C   
C      21 21      21 21      21 21
C    G(IK,JL) = G(IK,JL) + G(KI,JL) 
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
      NNM1O2(I)   =(I*(I-1))/2
C
C COMPUTE OFFSETS TO FIRST ELEMENT OF ALL G(XX,YY) IRREP POSITIONS
C  AND PUT IN FIRST CONTRIBUTION.
C
      LISTIJKL=110+ISPIN
      IOFFR=1
      DO 1000 IRREPR=1,NIRREP
       IOFFSR(IRREPR)=IOFFR
       IOFFSL(IRREPR)=IOFFR
       IOFFR=IOFFR+NNM1O2(POP(IRREPR,ISPIN))
1000  CONTINUE
      DO 100 IRREPDO=2,NIRREP
       DO 110 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPDO)
        IF(IRREPL.LT.IRREPR)THEN
         ILOGREC=IOFFSR(IRREPR)
         IFIRST =IOFFSL(IRREPL)
         NUML=POP(IRREPR,ISPIN)
         NUMK=POP(IRREPR,ISPIN)
         NUMJ=POP(IRREPL,ISPIN)
         NUMI=POP(IRREPL,ISPIN)
         MAXL=MOPOP(IRREPL)*MOPOP(IRREPR)
         MAXA=MOPOP(IRREPR)
         DO 120 INDXL=2,NUML
          DO 130 INDXK=1,INDXL-1
           INDL=INDXL
           INDK=INDXK
           CALL GETLST(BUF,ILOGREC,1,1,1,LISTIJKL)
           ILOGREC=ILOGREC+1
           IPOS=IFIRST
           DO 140 INDXJ=2,NUMJ
*VOCL LOOP,NOVREC
CDIR$ IVDEP
            DO 150 INDXI=1,INDXJ-1
             X=BUF(IPOS)*FIJKL
             IPOS=IPOS+1
             INDI=INDXI
             INDJ=INDXJ
             GINDXL=INDXF(INDK,INDI,MAXA)
             GINDXR=INDXF(INDL,INDJ,MAXA)
             GINDX0=INDXF(GINDXL,GINDXR,MAXL)
             GINDX =GINDX0+GOFF(IRREPR,IRREPL)
             GAMMA(GINDX)=GAMMA(GINDX)+X
             GINDX0=INDXF(GINDXR,GINDXL,MAXL)
             GINDX =GINDX0+GOFF(IRREPR,IRREPL)
             GAMMA(GINDX)=GAMMA(GINDX)+X
             GINDXL=INDXF(INDK,INDJ,MAXA)
             GINDXR=INDXF(INDL,INDI,MAXA)
             GINDX0=INDXF(GINDXL,GINDXR,MAXL)
             GINDX =GINDX0+GOFF(IRREPR,IRREPL)
             GAMMA(GINDX)=GAMMA(GINDX)-X
             GINDX0=INDXF(GINDXR,GINDXL,MAXL)
             GINDX =GINDX0+GOFF(IRREPR,IRREPL)
             GAMMA(GINDX)=GAMMA(GINDX)-X
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
C       21 21        12 21        12 12
C     G(KI,JL) <-- g<IJ,KL> = - g<IJ,LK>
C
      DO 160 IRREPDO=2,NIRREP
C
C COMPUTE OFFSETS TO ALL AB AB POSITIONS
C
      IOFFR=1
      IOFFL=1
      DO 1200 IRREPR=1,NIRREP
       IRREPL=DIRPRD(IRREPR,IRREPDO)
       IF(IRREPL.GT.IRREPR)GOTO 1200
       IOFFSR(IRREPR)=IOFFR
       IOFFSL(IRREPR)=IOFFR
       IOFFR=IOFFR+POP(IRREPL,ISPIN)*POP(IRREPR,ISPIN)
1200  CONTINUE
      DO 170 IRREPB=1,NIRREP
       IRREPA=DIRPRD(IRREPB,IRREPDO)
       IF(IRREPA.LT.IRREPB)THEN
        ILOGREC=IOFFSR(IRREPB)
        IFIRST=IOFFSL(IRREPB)
        NUML=POP(IRREPB,ISPIN)
        NUMK=POP(IRREPA,ISPIN)
        NUMJ=POP(IRREPB,ISPIN)
        NUMI=POP(IRREPA,ISPIN)
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
             GINDXL=INDXF(INDL,INDI,MAXB)
             GINDXR=INDXF(INDJ,INDK,MAXB)
             GINDX0=INDXF(GINDXL,GINDXR,MAXL)
             GINDX=GINDX0+GOFF(IRREPB,IRREPA)
             GAMMA(GINDX)=GAMMA(GINDX)-X
195         CONTINUE
190        CONTINUE
185      CONTINUE
180     CONTINUE
       ENDIF
170   CONTINUE
160   CONTINUE
      RETURN
      END