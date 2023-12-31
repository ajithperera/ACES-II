      SUBROUTINE G4AB42(GAMMA,BUF,GOFF,DOIRR,LENGAM,FAIBJ)
C
C THIS ROUTINE PROCESSES THE G(BI,aj) GAMMA LISTS AND FILLS
C  OUT THE ABCD SYMMETRY CASE MULLIKEN ORDERED GAMMA.
C
C THESE LISTS ARE A BIT TRICKY SINCE THEY ARE ALREADY STORED
C  IN MULLIKEN ORDER AND CONSEQUENTLY THE USUAL LOGIC DOES
C  NOT APPLY.
C
C WE MUST TAKE CARE OF THE VOVO, OVVO, OVOV AND VOOV PIECES
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION GAMMA(LENGAM),BUF(1),X,X1,X2 
      DOUBLE PRECISION FABIJ,FAIBJ,FABCD,FIJKL,FABCI,FIJKA
      INTEGER GOFF(8,8,8,8),IDID(8),IOFFSL(8)
      INTEGER IOFFSR(8),DOIRR(42,4)
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
C CALCULATE HOW MANY UNIQUE ABCD SYMMETRY SETS THERE ARE
C
      IF(NIRREP.EQ.8)THEN
       NUMABCD=42
      ELSEIF(NIRREP.EQ.4)THEN
       NUMABCD=3
      ENDIF
C     
      LISTBIAJ=118
      DO 100 ICOUNT=1,NUMABCD
       IRREP1M=DOIRR(ICOUNT,1)
       IRREP2M=DOIRR(ICOUNT,2)
       IRREP3M=DOIRR(ICOUNT,3)
       IRREP4M=DOIRR(ICOUNT,4)
       IOFFSET=GOFF(IRREP1M,IRREP2M,IRREP3M,IRREP4M)
       IRREPGET=DIRPRD(IRREP3M,IRREP4M)
C
       IOFFR=1
       IOFFL=1
       DO 2000 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREPGET)
        IOFFSR(IRREPA)=IOFFR
        IOFFSL(IRREPA)=IOFFL
        IOFFR=IOFFR+VRT(IRREPB,2)*POP(IRREPA,2)
        IOFFL=IOFFL+VRT(IRREPB,1)*POP(IRREPA,1)
2000   CONTINUE
C
C FIRST DO THE VOVO AND OVVO PARTS
C
       ILOGREC=IOFFSR(IRREP4M)
       IFIRST=IOFFSL(IRREP2M)
       NUMJ=POP(IRREP4M,2)
       NUMB=VRT(IRREP3M,2)
       NUMI2=POP(IRREP2M,1)
       NUMA1=VRT(IRREP1M,1)
       NUMI1=POP(IRREP1M,1)
       NUMA2=VRT(IRREP2M,1)
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
C
C VOVO PART
C
       DO 120 INDXJ=1,NUMJ
        DO 130 INDXB=1,NUMB
         INDJ=INDXJ
         INDB=INDXB+POP(IRREP3M,2)
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTBIAJ)
         ILOGREC=ILOGREC+1
         IPOS2=IOFFSL(IRREP2M)
         IPOS1=IOFFSL(IRREP1M)
         DO 140 INDXI=1,NUMI2
          DO 150 INDXA=1,NUMA1
           X=BUF(IPOS2)*FAIBJ
           IPOS2=IPOS2+1
           INDA=INDXA+POP(IRREP1M,1)
           INDI=INDXI
           GINDXL=INDXF(INDA,INDI,MOPOP(IRREP1M))
           GINDXR=INDXF(INDB,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
150       CONTINUE
140      CONTINUE
C
C OVVO PART
C
         DO 240 INDXI=1,NUMI1
          DO 250 INDXA=1,NUMA2
           X=BUF(IPOS1)*FAIBJ
           IPOS1=IPOS1+1
           INDA=INDXA+POP(IRREP2M,1)
           INDI=INDXI
           GINDXL=INDXF(INDI,INDA,MOPOP(IRREP1M))
           GINDXR=INDXF(INDB,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
250       CONTINUE
240      CONTINUE
130     CONTINUE
120    CONTINUE
C
C NOW FINISH WITH THE OVOV AND VOOV PARTS
C
       ILOGREC=IOFFSR(IRREP3M)
       IFIRST=IOFFSL(IRREP1M)
       NUMJ=POP(IRREP3M,2)
       NUMB=VRT(IRREP4M,2)
       NUMI1=POP(IRREP1M,1)
       NUMA2=VRT(IRREP2M,1)
       NUMI2=POP(IRREP2M,1)
       NUMA1=VRT(IRREP1M,1)
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
C
C OVOV PART
C
       DO 320 INDXJ=1,NUMJ
        DO 330 INDXB=1,NUMB
         INDJ=INDXJ
         INDB=INDXB+POP(IRREP4M,2)
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTBIAJ)
         ILOGREC=ILOGREC+1
         IPOS2=IOFFSL(IRREP2M)
         IPOS1=IOFFSL(IRREP1M)
         DO 340 INDXI=1,NUMI1
          DO 350 INDXA=1,NUMA2
           X=BUF(IPOS1)*FAIBJ
           IPOS1=IPOS1+1
           INDA=INDXA+POP(IRREP2M,1)
           INDI=INDXI
           GINDXL=INDXF(INDI,INDA,MOPOP(IRREP1M))
           GINDXR=INDXF(INDJ,INDB,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
350       CONTINUE
340      CONTINUE
C
C VOOV PART
C
         DO 360 INDXI=1,NUMI2
          DO 370 INDXA=1,NUMA1
           X=BUF(IPOS2)*FAIBJ
           IPOS2=IPOS2+1
           INDA=INDXA+POP(IRREP1M,1)
           INDI=INDXI
           GINDXL=INDXF(INDA,INDI,MOPOP(IRREP1M))
           GINDXR=INDXF(INDJ,INDB,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
370       CONTINUE
360      CONTINUE
330     CONTINUE
320    CONTINUE
C
100   CONTINUE
      RETURN
      END
