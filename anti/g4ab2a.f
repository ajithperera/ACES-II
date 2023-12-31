      SUBROUTINE G4AB2A(GAMMA,BUF,GOFF,DOIRR,LENGAM,FIJKA)
C
C THIS ROUTINE PROCESSES THE G<Ij,Ka> AND G<Ij,Ak> GAMMA LISTS 
C  AND FILLS OUT THE ABCD SYMMETRY CASE MULLIKEN ORDERED GAMMA FOR 
C  SPIN CASE AABB.  THIS IS FOR UHF CALCULATIONS ONLY!
C
C THE COMMENTS IN THIS ROUTINE REFER TO THE <Ab,Ci> CASE, BUT
C  ARE CONSISTENT WITH THE CODE PROVIDED THE FOLLOWING SUBSTITUTIONS
C  ARE MADE:
C
C             I<-A J<-B K<-C A<-I
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
C
C CALCULATE HOW MANY UNIQUE ABCD SYMMETRY SETS THERE ARE
C
      IF(NIRREP.EQ.8)THEN
       NUMABCD=42
      ELSEIF(NIRREP.EQ.4)THEN
       NUMABCD=3
      ENDIF
C     
      DO 100 ICOUNT=1,NUMABCD
       IRREP1M=DOIRR(ICOUNT,1)
       IRREP2M=DOIRR(ICOUNT,2)
       IRREP3M=DOIRR(ICOUNT,3)
       IRREP4M=DOIRR(ICOUNT,4)
       IOFFSET=GOFF(IRREP1M,IRREP2M,IRREP3M,IRREP4M)
C
C FIRST DO LIST 110 
C
C          12 34        13 24
C        G(AC,bi) <-- g<Ab,Ci>
C
C          12 34        24 13
C        G(CA,ib) <-- g<Ab,Ci>
C
C          12 34        23 14
C        G(CA,bi) <-- g<Ab,Ci>
C
C          12 34        14 23
C        G(AC,ib) <-- g<Ab,Ci>
C
       LISTIJKA=110
       IRREPGET=DIRPRD(IRREP2M,IRREP4M)
C
       IOFFR=1
       IOFFL=1
       DO 2000 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREPGET)
        IOFFSR(IRREPA)=IOFFR
        IOFFSL(IRREPA)=IOFFL
        IOFFR=IOFFR+POP(IRREPB,1)*VRT(IRREPA,2)
        IOFFL=IOFFL+POP(IRREPB,1)*POP(IRREPA,2)
2000   CONTINUE
C
C          12 34        13 24
C        G(AC,bi) <-- g<Ab,Ci>
C
       ILOGREC=IOFFSR(IRREP4M)
       IFIRST=IOFFSL(IRREP3M)
       NUMA=VRT(IRREP4M,2)
       NUMK=POP(IRREP2M,1)
       NUMJ=POP(IRREP3M,2)
       NUMI=POP(IRREP1M,1) 
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
       DO 120 INDXA=1,NUMA
        DO 130 INDXK=1,NUMK
         INDA=INDXA+POP(IRREP4M,2)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP3M)
         DO 140 INDXJ=1,NUMJ
          DO 150 INDXI=1,NUMI
           X=BUF(IPOS)*FIJKA
           IPOS=IPOS+1
           INDJ=INDXJ
           INDI=INDXI
           GINDXL=INDXF(INDI,INDK,MOPOP(IRREP1M))
           GINDXR=INDXF(INDJ,INDA,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
150       CONTINUE
140      CONTINUE
130     CONTINUE
120    CONTINUE
C
C          12 34        24 13
C        G(CA,ib) <-- g<Ab,Ci>
C
       ILOGREC=IOFFSR(IRREP3M)
       IFIRST=IOFFSL(IRREP4M)
       NUMA=VRT(IRREP3M,2)
       NUMK=POP(IRREP1M,1) 
       NUMJ=POP(IRREP4M,2)
       NUMI=POP(IRREP2M,1)
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
       DO 420 INDXA=1,NUMA
        DO 430 INDXK=1,NUMK
         INDA=INDXA+POP(IRREP3M,2)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP4M)
         DO 440 INDXJ=1,NUMJ
          DO 450 INDXI=1,NUMI
           X=BUF(IPOS)*FIJKA
           IPOS=IPOS+1
           INDJ=INDXJ
           INDI=INDXI
           GINDXL=INDXF(INDK,INDI,MOPOP(IRREP1M))
           GINDXR=INDXF(INDA,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
450       CONTINUE
440      CONTINUE
430     CONTINUE
420    CONTINUE
C
C          12 34        23 14
C        G(CA,bi) <-- g<Ab,Ci>
C
       IRREPGET=DIRPRD(IRREP1M,IRREP4M)
       IOFFR=1
       IOFFL=1
C
       DO 1000 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREPGET)
        IOFFSR(IRREPA)=IOFFR
        IOFFSL(IRREPA)=IOFFL
        IOFFR=IOFFR+POP(IRREPB,1)*VRT(IRREPA,2)
        IOFFL=IOFFL+POP(IRREPB,1)*POP(IRREPA,2)
1000   CONTINUE
C
       ILOGREC=IOFFSR(IRREP4M)
       IFIRST=IOFFSL(IRREP3M)
       NUMA=VRT(IRREP4M,2)
       NUMK=POP(IRREP1M,1) 
       NUMJ=POP(IRREP3M,2)
       NUMI=POP(IRREP2M,1)
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
       DO 220 INDXA=1,NUMA
        DO 230 INDXK=1,NUMK
         INDA=INDXA+POP(IRREP4M,2)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP3M)
         DO 240 INDXJ=1,NUMJ
          DO 250 INDXI=1,NUMI
           X=BUF(IPOS)*FIJKA
           IPOS=IPOS+1
           INDJ=INDXJ
           INDI=INDXI
           GINDXL=INDXF(INDK,INDI,MOPOP(IRREP1M))
           GINDXR=INDXF(INDJ,INDA,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
250       CONTINUE
240      CONTINUE
230     CONTINUE
220    CONTINUE
C
C          12 34        14 23
C        G(IK,aj) <-- g<Ij,Ka>
C
       ILOGREC=IOFFSR(IRREP3M)
       IFIRST=IOFFSL(IRREP4M)
       NUMA=VRT(IRREP3M,2)
       NUMK=POP(IRREP2M,1) 
       NUMJ=POP(IRREP4M,2)
       NUMI=POP(IRREP1M,1)
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
       DO 320 INDXA=1,NUMA
        DO 330 INDXK=1,NUMK
         INDA=INDXA+POP(IRREP3M,2)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP4M)
         DO 340 INDXJ=1,NUMJ
          DO 350 INDXI=1,NUMI
           X=BUF(IPOS)*FIJKA
           IPOS=IPOS+1
           INDJ=INDXJ
           INDI=INDXI
           GINDXL=INDXF(INDI,INDK,MOPOP(IRREP1M))
           GINDXR=INDXF(INDA,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
350       CONTINUE
340      CONTINUE
330     CONTINUE
320    CONTINUE
C
C NOW DO LIST 109
C
C          12 34        13 24
C        G(AI,bc) <-- g<Ab,Ic>
C
C          12 34        24 13
C        G(IA,cb) <-- g<Ab,Ic>
C
C          12 34        23 14
C        G(IA,bc) <-- g<Ab,Ic>
C
C          12 34        14 23
C        G(AI,cb) <-- g<Ab,Ic>
C
       LISTIJKA=109
       IRREPGET=DIRPRD(IRREP2M,IRREP4M)
C
       IOFFR=1
       IOFFL=1
       DO 4000 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREPGET)
        IOFFSR(IRREPA)=IOFFR
        IOFFSL(IRREPA)=IOFFL
        IOFFR=IOFFR+VRT(IRREPB,1)*POP(IRREPA,2)
        IOFFL=IOFFL+POP(IRREPB,1)*POP(IRREPA,2)
4000   CONTINUE
C
C          12 34        13 24
C        G(AI,bc) <-- g<Ab,Ic>
C
       ILOGREC=IOFFSR(IRREP4M)
       IFIRST=IOFFSL(IRREP3M)
       NUMK=POP(IRREP4M,2)
       NUMA=VRT(IRREP2M,1)
       NUMJ=POP(IRREP3M,2)
       NUMI=POP(IRREP1M,1) 
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
       DO 1120 INDXK=1,NUMK
        DO 1130 INDXA=1,NUMA
         INDA=INDXA+POP(IRREP2M,1)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP3M)
         DO 1140 INDXJ=1,NUMJ
          DO 1150 INDXI=1,NUMI
           X=BUF(IPOS)*FIJKA
           IPOS=IPOS+1
           INDJ=INDXJ
           INDI=INDXI
           GINDXL=INDXF(INDI,INDA,MOPOP(IRREP1M))
           GINDXR=INDXF(INDJ,INDK,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
1150       CONTINUE
1140      CONTINUE
1130     CONTINUE
1120    CONTINUE
C
C          12 34        24 13
C        G(IA,cb) <-- g<Ab,Ic>
C
       ILOGREC=IOFFSR(IRREP3M)
       IFIRST=IOFFSL(IRREP4M)
       NUMK=POP(IRREP3M,2) 
       NUMA=VRT(IRREP1M,1)
       NUMJ=POP(IRREP4M,2)
       NUMI=POP(IRREP2M,1)
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
       DO 1420 INDXK=1,NUMK
        DO 1430 INDXA=1,NUMA
         INDA=INDXA+POP(IRREP1M,1)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP4M)
         DO 1440 INDXJ=1,NUMJ
          DO 1450 INDXI=1,NUMI
           X=BUF(IPOS)*FIJKA
           IPOS=IPOS+1
           INDJ=INDXJ
           INDI=INDXI
           GINDXL=INDXF(INDA,INDI,MOPOP(IRREP1M))
           GINDXR=INDXF(INDK,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
1450       CONTINUE
1440      CONTINUE
1430     CONTINUE
1420    CONTINUE
C
C          12 34        23 14
C        G(IA,bc) <-- g<Ab,Ic>
C
       IRREPGET=DIRPRD(IRREP1M,IRREP4M)
       IOFFR=1
       IOFFL=1
C
       DO 5000 IRREPA=1,NIRREP
        IRREPB=DIRPRD(IRREPA,IRREPGET)
        IOFFSR(IRREPA)=IOFFR
        IOFFSL(IRREPA)=IOFFL
        IOFFR=IOFFR+VRT(IRREPB,1)*POP(IRREPA,2)
        IOFFL=IOFFL+POP(IRREPB,1)*POP(IRREPA,2)
5000   CONTINUE
C
       ILOGREC=IOFFSR(IRREP4M)
       IFIRST=IOFFSL(IRREP3M)
       NUMK=POP(IRREP4M,2) 
       NUMA=VRT(IRREP1M,1)
       NUMJ=POP(IRREP3M,2)
       NUMI=POP(IRREP2M,1)
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
       DO 1220 INDXK=1,NUMK
        DO 1230 INDXA=1,NUMA
         INDA=INDXA+POP(IRREP1M,1)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP3M)
         DO 1240 INDXJ=1,NUMJ
          DO 1250 INDXI=1,NUMI
           X=BUF(IPOS)*FIJKA
           IPOS=IPOS+1
           INDJ=INDXJ
           INDI=INDXI
           GINDXL=INDXF(INDA,INDI,MOPOP(IRREP1M))
           GINDXR=INDXF(INDJ,INDK,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
1250       CONTINUE
1240      CONTINUE
1230     CONTINUE
1220    CONTINUE
C
C          12 34        14 23
C        G(AI,cb) <-- g<Ab,Ic>
C
       ILOGREC=IOFFSR(IRREP3M)
       IFIRST=IOFFSL(IRREP4M)
       NUMK=POP(IRREP3M,2) 
       NUMA=VRT(IRREP2M,1)
       NUMJ=POP(IRREP4M,2)
       NUMI=POP(IRREP1M,1)
       MAXL=MOPOP(IRREP1M)*MOPOP(IRREP2M)
       DO 1320 INDXK=1,NUMK
        DO 1330 INDXA=1,NUMA
         INDA=INDXA+POP(IRREP2M,1)
         INDK=INDXK
         CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTIJKA)
         ILOGREC=ILOGREC+1
         IPOS=IOFFSL(IRREP4M)
         DO 1340 INDXJ=1,NUMJ
          DO 1350 INDXI=1,NUMI
           X=BUF(IPOS)*FIJKA
           IPOS=IPOS+1
           INDJ=INDXJ
           INDI=INDXI
           GINDXL=INDXF(INDI,INDA,MOPOP(IRREP1M))
           GINDXR=INDXF(INDK,INDJ,MOPOP(IRREP3M))
           GINDX0=INDXF(GINDXL,GINDXR,MAXL)
           GINDX =GINDX0+IOFFSET
           GAMMA(GINDX)=GAMMA(GINDX)+X
1350       CONTINUE
1340      CONTINUE
1330     CONTINUE
1320    CONTINUE
C
100   CONTINUE
      RETURN
      END
