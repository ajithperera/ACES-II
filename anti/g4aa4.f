      SUBROUTINE G4AA4(GAMMA,BUF,GOFF,DOIRR,LENGAM,ISPIN,FAIBJ)
C
C THIS ROUTINE PROCESSES THE G<AI,BJ> GAMMA LISTS AND FILLS
C  OUT THE ABCD SYMMETRY CASE MULLIKEN ORDERED GAMMA FOR 
C  SPIN CASE AAAA.  THIS IS FOR UHF CALCULATIONS ONLY!
C
C          12 34        13 24
C        G(AB,IJ) <-- g<AI,BJ>
C
C          12 34        31 42
C        G(IJ,AB) <-- g<AI,BJ>
C
C          12 34        23 14
C        G(BA,IJ) <-- g<AI,BJ>
C
C          12 34        32 41
C        G(JI,AB) <-- g<AI,BJ>
C
C          12 34        13 24        13 42
C        G(AJ,IB) <-- g<AI,JB> = - g<AI,BJ>
C
C          12 34        13 24        31 24
C        G(IA,BJ) <-- g<IB,AJ> = - g<BI,AJ>
C
C          12 34        14 23        14 32        32 14
C        G(BI,AJ) <-- g<BJ,IA> = - g<BJ,AI> = - g<AI,BJ>
C
C          12 34        14 23        41 23        23 41
C        G(JA,IB) <-- g<JB,AI> = - g<BJ,AI> = - g<AI,BJ>
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION GAMMA(LENGAM),BUF(1),X,X1,X2 
      DOUBLE PRECISION FABIJ,FIJKL,FABCD,FAIBJ,FABCI,FIJKA
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
      NNM1O2(I)   =(I*(I-1))/2
C
C CALCULATE HOW MANY UNIQUE ABCD SYMMETRY SETS THERE ARE
C
      IF(NIRREP.EQ.8)THEN
       NUMABCD=42
      ELSEIF(NIRREP.EQ.4)THEN
       NUMABCD=3
      ENDIF
C     
      LISTAIBJ=122+ISPIN
      DO 100 ICOUNT=1,NUMABCD
       IRREP1M=DOIRR(ICOUNT,1)
       IRREP2M=DOIRR(ICOUNT,2)
       IRREP3M=DOIRR(ICOUNT,3)
       IRREP4M=DOIRR(ICOUNT,4)
C
       IOFFSET=GOFF(IRREP1M,IRREP2M,IRREP3M,IRREP4M)
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
       MAXL1=IRREP1M
       MAXL2=IRREP3M
C
C         12 34        13 24
C       G(AB,IJ) <-- g<AI,BJ>    
C
       IRREPGET=DIRPRD(IRREP2M,IRREP4M)
C
       IOFFR=1
       IOFFL=1
       DO 2000 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREPGET)
        IOFFSR(IRREPA)=IOFFR
        IOFFSL(IRREPA)=IOFFL
        IOFFR=IOFFR+VRT(IRREPB,ISPIN)*POP(IRREPA,ISPIN)
        IOFFL=IOFFL+VRT(IRREPB,ISPIN)*POP(IRREPA,ISPIN)
2000   CONTINUE
C
       ILOGREC=IOFFSR(IRREP4M)
       IFIRST=IOFFSL(IRREP3M)
       NUMJ=POP(IRREP4M,ISPIN)
       NUMB=VRT(IRREP2M,ISPIN)
       NUMI=POP(IRREP3M,ISPIN)
       NUMA=VRT(IRREP1M,ISPIN) 
       NUMI2=POP(IRREP1M,ISPIN)
       NUMA2=VRT(IRREP3M,ISPIN) 
       DO 120 INDXJ=1,NUMJ
        DO 130 INDXB=1,NUMB
         INDJ=INDXJ 
         INDB=INDXB+POP(IRREP2M,ISPIN)
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTAIBJ)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP3M)
         DO 140 INDXI=1,NUMI
          DO 150 INDXA=1,NUMA
           X=BUF(IPOS)*FAIBJ
           IPOS=IPOS+1
           INDI=INDXI
           INDA=INDXA+POP(IRREP1M,ISPIN)
           GINDXL=INDXF(INDA,INDB,MOPOP(IRREP1M))
           GINDXR=INDXF(INDI,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
150       CONTINUE
140      CONTINUE
C
C          12 34        13 24        31 24
C        G(IB,AJ) <-- g<IA,BJ> = - g<AI,BJ>
C
         IPOS=IOFFSL(IRREP1M)
         DO 141 INDXI=1,NUMI2
          DO 151 INDXA=1,NUMA2
           X=BUF(IPOS)*FAIBJ
           IPOS=IPOS+1
           INDI=INDXI
           INDA=INDXA+POP(IRREP3M,ISPIN)
           GINDXL=INDXF(INDI,INDB,MOPOP(IRREP1M))
           GINDXR=INDXF(INDA,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)-X
151       CONTINUE
141      CONTINUE
130     CONTINUE
120    CONTINUE
C
C         12 34        31 42
C       G(IJ,AB) <-- g<AI,BJ>  
C
       ILOGREC=IOFFSR(IRREP2M)
       IFIRST=IOFFSL(IRREP1M)
       NUMJ=POP(IRREP2M,ISPIN)
       NUMB=VRT(IRREP4M,ISPIN)
       NUMI=POP(IRREP1M,ISPIN)
       NUMA=VRT(IRREP3M,ISPIN) 
       NUMI2=POP(IRREP3M,ISPIN)
       NUMA2=VRT(IRREP1M,ISPIN) 
       DO 220 INDXJ=1,NUMJ
        DO 230 INDXB=1,NUMB
         INDJ=INDXJ 
         INDB=INDXB+POP(IRREP4M,ISPIN)
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTAIBJ)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP1M)
         DO 240 INDXI=1,NUMI
          DO 250 INDXA=1,NUMA
           X=BUF(IPOS)*FAIBJ
           IPOS=IPOS+1
           INDI=INDXI
           INDA=INDXA+POP(IRREP3M,ISPIN)
           GINDXL=INDXF(INDI,INDJ,MOPOP(IRREP1M))
           GINDXR=INDXF(INDA,INDB,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
250       CONTINUE
240      CONTINUE
C
C          12 34        13 24        13 42
C        G(AJ,IB) <-- g<AI,JB> = - g<AI,BJ>
C
         IPOS=IOFFSL(IRREP3M)
         DO 241 INDXI=1,NUMI2
          DO 251 INDXA=1,NUMA2
           X=BUF(IPOS)*FAIBJ
           IPOS=IPOS+1
           INDI=INDXI
           INDA=INDXA+POP(IRREP1M,ISPIN)
           GINDXL=INDXF(INDA,INDJ,MOPOP(IRREP1M))
           GINDXR=INDXF(INDI,INDB,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)-X
251       CONTINUE
241      CONTINUE
C
230     CONTINUE
220    CONTINUE
C
C         12 34        23 14
C       G(BA,IJ) <-- g<AI,BJ>
C
       IRREPGET=DIRPRD(IRREP1M,IRREP4M)
C
       IOFFR=1
       IOFFL=1
       DO 3000 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREPGET)
        IOFFSR(IRREPA)=IOFFR
        IOFFSL(IRREPA)=IOFFL
        IOFFR=IOFFR+VRT(IRREPB,ISPIN)*POP(IRREPA,ISPIN)
        IOFFL=IOFFL+VRT(IRREPB,ISPIN)*POP(IRREPA,ISPIN)
3000   CONTINUE
C
       ILOGREC=IOFFSR(IRREP4M)
       IFIRST=IOFFSL(IRREP3M)
       NUMJ=POP(IRREP4M,ISPIN)
       NUMB=VRT(IRREP1M,ISPIN)
       NUMI=POP(IRREP3M,ISPIN)
       NUMA=VRT(IRREP2M,ISPIN) 
       NUMI2=POP(IRREP2M,ISPIN)
       NUMA2=VRT(IRREP3M,ISPIN) 
       DO 320 INDXJ=1,NUMJ
        DO 330 INDXB=1,NUMB
         INDJ=INDXJ 
         INDB=INDXB+POP(IRREP1M,ISPIN)
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTAIBJ)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP3M)
         DO 340 INDXI=1,NUMI
          DO 350 INDXA=1,NUMA
           X=BUF(IPOS)*FAIBJ
           IPOS=IPOS+1
           INDI=INDXI
           INDA=INDXA+POP(IRREP2M,ISPIN)
           GINDXL=INDXF(INDB,INDA,MOPOP(IRREP1M))
           GINDXR=INDXF(INDI,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
350       CONTINUE
340      CONTINUE
C
C          12 34        14 23        14 32        32 14
C        G(BI,AJ) <-- g<BJ,IA> = - g<BJ,AI> = - g<AI,BJ>
C
         IPOS=IOFFSL(IRREP2M)
         DO 341 INDXI=1,NUMI2
          DO 351 INDXA=1,NUMA2
           X=BUF(IPOS)*FAIBJ
           IPOS=IPOS+1
           INDI=INDXI
           INDA=INDXA+POP(IRREP3M,ISPIN)
           GINDXL=INDXF(INDB,INDI,MOPOP(IRREP1M))
           GINDXR=INDXF(INDA,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)-X
351       CONTINUE
341      CONTINUE
330     CONTINUE
320    CONTINUE
C
C          12 34        32 41
C        G(JI,AB) <-- g<AI,BJ>
C
       ILOGREC=IOFFSR(IRREP1M)
       IFIRST=IOFFSL(IRREP2M)
       NUMJ=POP(IRREP1M,ISPIN)
       NUMB=VRT(IRREP4M,ISPIN)
       NUMI=POP(IRREP2M,ISPIN)
       NUMA=VRT(IRREP3M,ISPIN) 
       NUMI2=POP(IRREP3M,ISPIN)
       NUMA2=VRT(IRREP2M,ISPIN) 
       DO 420 INDXJ=1,NUMJ
        DO 430 INDXB=1,NUMB
         INDJ=INDXJ 
         INDB=INDXB+POP(IRREP4M,ISPIN)
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTAIBJ)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP2M)
         DO 440 INDXI=1,NUMI
          DO 450 INDXA=1,NUMA
           X=BUF(IPOS)*FAIBJ
           IPOS=IPOS+1
           INDI=INDXI
           INDA=INDXA+POP(IRREP3M,ISPIN)
           GINDXL=INDXF(INDJ,INDI,MOPOP(IRREP1M))
           GINDXR=INDXF(INDA,INDB,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
450       CONTINUE
440      CONTINUE
C
C          12 34        14 23        41 23        23 41
C        G(JA,IB) <-- g<JB,AI> = - g<BJ,AI> = - g<AI,BJ>
C
         IPOS=IOFFSL(IRREP3M)
         DO 441 INDXI=1,NUMI2
          DO 451 INDXA=1,NUMA2
           X=BUF(IPOS)*FAIBJ
           IPOS=IPOS+1
           INDI=INDXI
           INDA=INDXA+POP(IRREP2M,ISPIN)
           GINDXL=INDXF(INDJ,INDA,MOPOP(IRREP1M))
           GINDXR=INDXF(INDI,INDB,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)-X
451       CONTINUE
441      CONTINUE
430     CONTINUE
420    CONTINUE
C   
100   CONTINUE
      RETURN
      END
