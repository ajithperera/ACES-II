      SUBROUTINE G2AB1A(GAMMA,BUF,GOFF,LENGAM,FIJKL)
C
C THIS ROUTINE PROCESSES THE G(Ij,Kl) GAMMAS AND AUGMENTS THE 
C  AABB MULLIKEN ORDERED MO GAMMA LIST FOR SPIN CASE AABB.  THIS
C  CODE USED IN UHF CALCULATIONS.  THE FOLLOWING SYMMETRIZATION
C  IS PERFORMED AS WELL.
C      
C           11 22        12 12      12 12
C         G(IK,jl) <-- g<Ij,Kl> + g<Kj,Il>
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION GAMMA(LENGAM),BUF(1),X,X1,X2 
      DOUBLE PRECISION FABIJ,FAIBJ,FABCD,FIJKL,FABCI,FIJKA
      INTEGER GOFF(8,8),IOFFSR(8),IOFFSL(8)
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
      LISTIJKL=113
      DO 100 IRREPGET=2,NIRREP
       IOFFR=1
       IOFFL=1
       DO 1100 IRREP1=1,NIRREP
        IRREP2=DIRPRD(IRREP1,IRREPGET)
        IOFFSR(IRREP1)=IOFFR
        IOFFSL(IRREP1)=IOFFL
        IOFFR=IOFFR+POP(IRREP2,1)*POP(IRREP1,2)
        IOFFL=IOFFL+POP(IRREP2,1)*POP(IRREP1,2)
1100   CONTINUE
       DO 110 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPGET)
        IF(IRREPL.GT.IRREPR)GOTO 110
        ILOGREC=IOFFSR(IRREPR)
        IFIRST =IOFFSL(IRREPR)-1
        NUML=POP(IRREPR,2)
        NUMK=POP(IRREPL,1)
        NUMJ=POP(IRREPR,2)
        NUMI=POP(IRREPL,1)
        MAXL=NNP1O2(MOPOP(IRREPL))
        DO 120 INDXL=1,NUML
         DO 130 INDXK=1,NUMK
          INDL=INDXL
          INDK=INDXK
          CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTIJKL)
          ILOGREC=ILOGREC+1
          IPOS=IFIRST
          DO 140 INDXJ=1,INDXL
           INDJ=INDXJ
           GINDXR=INDXT(INDJ,INDL)
           DO 150 INDXI=1,INDXK
            IPOS=IFIRST+INDXF(INDXI,INDXJ,NUMI)
            X=BUF(IPOS)*FIJKL
            INDI=INDXI
            GINDXL=INDXT(INDI,INDK)
            GINDX0=INDXF(GINDXL,GINDXR,MAXL)
            GINDX =GINDX0+GOFF(IRREPL,IRREPR)
            GAMMA(GINDX)=GAMMA(GINDX)+X
150        CONTINUE
           DO 151 INDXI=INDXK,NUMI
            IPOS=IFIRST+INDXF(INDXI,INDXJ,NUMI)
            X=BUF(IPOS)*FIJKL
            INDI=INDXI
            GINDXL=INDXT(INDK,INDI)
            GINDX0=INDXF(GINDXL,GINDXR,MAXL)
            GINDX=GINDX0+GOFF(IRREPL,IRREPR)
            GAMMA(GINDX)=GAMMA(GINDX)+X
151        CONTINUE
140       CONTINUE
130      CONTINUE
120     CONTINUE
110    CONTINUE
100   CONTINUE
      RETURN
      END
