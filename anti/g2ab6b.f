      SUBROUTINE G2AB6B(GAMMA,BUF,BUF1,BUF2,GOFF,LENGAM,FABCD)
C
C THIS ROUTINE PROCESSES THE G(Ab,Cd) GAMMAS AND AUGMENTS THE 
C  AABB MULLIKEN ORDERED MO GAMMA LIST FOR SPIN CASE BBAA.
C
C            11 22        21 21
C          G(bd,AC) <-- g<Ab,Cd>
C
C            11 22        21 21
C          G(db,AC) <-- g<Ab,Cd>
C
C IN ADDITION, THE FOLLOWING SYMMETRIZATION IS REQUIRED
C
C      G(ab,CD) = G(ab,CD) + G(ba,CD)
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION GAMMA(LENGAM),BUF(*),X,X1,X2 
      DOUBLE PRECISION FABIJ,FAIBJ,FABCD,FIJKL,FABCI,FIJKA
      DOUBLE PRECISION BUF1(*),BUF2(*),ONE,TWO
      INTEGER GOFF(8,8),IOFFSR(8),IOFFSL(8)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      LOGICAL MBPT4
      LOGICAL GABCD,DABCD
      LOGICAL TDA,EOM 
C
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /MOPOPS/ MOPOP(8)
      COMMON /IABCD/ IOFFABCD
      COMMON /ABCD/  GABCD,DABCD
      COMMON /EXCITE/TDA,EOM
      COMMON /STATSYM/IRREPX
C
      DATA ONE,TWO/1.D0,2.D0/
C
      INDXF(I,J,N)=I+(J-1)*N
      INDXT(I,J)  =I+(J*(J-1))/2
      NNP1O2(I)   =(I*(I+1))/2
C
      NPASS=1
      IF(DABCD.AND.EOM) NPASS=2
C
      LISTABCD=IOFFABCD+133
C
      DO 100 IRREPGET=2,NIRREP
C
       DO 101 IPASS=1,NPASS
C
       IF(DABCD) THEN
C
        IF(IPASS.EQ.1) THEN
C
         NUMDSZ=IRPDPD(IRREPGET,ISYTYP(2,146))
         DISDSZ=IRPDPD(IRREPGET,ISYTYP(1,146))
         CALL GETLST(BUF1,1,NUMDSZ,1,IRREPGET,156)
         CALL GETLST(BUF2,1,NUMDSZ,1,IRREPGET,146)
         IF(CCSD) THEN
          I0TA=1
          I0TB=1+NT(1)
          CALL GETLST(BUF(I0TA),1,1,1,1,157)
          CALL GETLST(BUF(I0TB),1,1,1,2,157)
          CALL FTAU(BUF1,BUF(I0TA),BUF(I0TB),DISDSZ,NUMDSZ,
     &              POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),
     &              IRREPGET,3,ONE)
         ENDIF
C
         MBPT4=M4DQ.OR.M4SDQ.OR.M4SDTQ 
         IF(MBPT4) THEN
          CALL SSCAL(NUMDSZ*DISDSZ,TWO,BUF2,1) 
          CALL SAXPY(NUMDSZ*DISDSZ,ONE,BUF1,1,BUF2,1)
         ENDIF
C
        ELSE
C
C SECOND PASS, ONLY FOR EOM-CCSD REQUIRED
C
         IRREPGET2=DIRPRD(IRREPGET,IRREPX)
         NUMDSZ=IRPDPD(IRREPGET2,ISYTYP(2,463))
         DISDSZ=IRPDPD(IRREPGET,ISYTYP(1,463))
C
         CALL GETLST(BUF1,1,NUMDSZ,1,IRREPGET2,463)
         CALL GETLST(BUF2,1,NUMDSZ,1,IRREPGET2,446)
C
         I0TA=1
         I0TB=1+NT(1)
         I0RA=I0TB+NT(2)
         I0RB=I0RA+IRPDPD(IRREPX,9)
         CALL GETLST(BUF(I0TA),1,1,1,1,157)
         CALL GETLST(BUF(I0RA),1,1,1,3,490)
         CALL GETLST(BUF(I0TB),1,1,1,2,157)
         CALL GETLST(BUF(I0RB),1,1,1,4,490)
C
         CALL DTAU(IRREPGET,IRREPGET2,1,IRREPX,BUF1,
     &             BUF(I0TA),BUF(I0TB),BUF(I0RA),BUF(I0RB),
     &             3,ONE)
C
        ENDIF
C
       ENDIF
C
       IOFFR=1
       IOFFL=1
       DO 1100 IRREP1=1,NIRREP
        IRREP2=DIRPRD(IRREP1,IRREPGET)
        IOFFSR(IRREP1)=IOFFR
        IOFFSL(IRREP1)=IOFFL
        IOFFR=IOFFR+VRT(IRREP2,1)*VRT(IRREP1,2)
        IOFFL=IOFFL+VRT(IRREP2,1)*VRT(IRREP1,2)
1100   CONTINUE
       DO 110 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPGET)
        IF(IRREPL.LT.IRREPR)GOTO 110
        ILOGREC=IOFFSR(IRREPR)
        IFIRST =IOFFSL(IRREPR)-1
        NUMD=VRT(IRREPR,2)
        NUMC=VRT(IRREPL,1)
        NUMB=VRT(IRREPR,2)
        NUMA=VRT(IRREPL,1)
        MAXL=NNP1O2(MOPOP(IRREPR))
        DO 120 INDXD=1,NUMD
         DO 130 INDXC=1,NUMC
          INDD=INDXD+POP(IRREPR,2)
          INDC=INDXC+POP(IRREPL,1)
          IF(.NOT.DABCD) THEN
C
C READ G(AB,CD) FROM DISK
C
           CALL GETLST(BUF,ILOGREC,1,1,IRREPGET,LISTABCD)
C
          ELSE
C
C CALCULATE G(AB,CD) ON THE FLY
C
           CALL DIRG2(BUF,BUF1,BUF2,BUF1,BUF2,NUMDSZ,DISDSZ,
     &                ILOGREC,IRREPGET,1)
C
          ENDIF
          ILOGREC=ILOGREC+1
          DO 140 INDXB=1,INDXD
           INDB=INDXB+POP(IRREPR,2)
           GINDXR=INDXT(INDB,INDD)
           DO 150 INDXA=1,INDXC
            IPOS=IFIRST+INDXF(INDXA,INDXB,NUMA)
            X=BUF(IPOS)*FABCD
            INDA=INDXA+POP(IRREPL,1)
            GINDXL=INDXT(INDA,INDC)
            GINDX0=INDXF(GINDXR,GINDXL,MAXL)
            GINDX =GINDX0+GOFF(IRREPR,IRREPL)
            GAMMA(GINDX)=GAMMA(GINDX)+X
150        CONTINUE
           DO 151 INDXA=INDXC,NUMA
            IPOS=IFIRST+INDXF(INDXA,INDXB,NUMA)
            X=BUF(IPOS)*FABCD
            INDA=INDXA+POP(IRREPL,1)
            GINDXL=INDXT(INDC,INDA)
            GINDX0=INDXF(GINDXR,GINDXL,MAXL)
            GINDX=GINDX0+GOFF(IRREPR,IRREPL)
            GAMMA(GINDX)=GAMMA(GINDX)+X
151        CONTINUE
140       CONTINUE
130      CONTINUE
120     CONTINUE
110    CONTINUE
101   CONTINUE
100   CONTINUE
      RETURN
      END
