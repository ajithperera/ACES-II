      SUBROUTINE OST2IIAI2(ICORE,MAXCOR,IUHF,W,LAMBDA,WL)
C
C
C     THIS ROUTINE CALCULATE INACTIVE,INACTIVE-ACTIVE,INACTIVE T2 CONTRIBUTION
C
C     T2(Ij,Am)=-W * TAU(iJ,mN)*T1(n,a)
C     T2(Ij,Na)=-W * TAU(iJ,mN)*T1(M,A)
C
CEND
C
CPROGRAMED PS NOV/92
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,POP,VRT,DISSYT2
      CHARACTER*4 STRING
      LOGICAL LAMBDA
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
      INTEGER POPA,POPI,VRTA,VRTI
      COMMON /FSSYM/ POPA(8,2),POPI(8,2),VRTA(8,2),VRTI(8,2),
     &               NTAN(2),NTNA(2),NTAA(2),
     &               NTIN(2),NTNI(2),NTII(2),
     &               NTAI(2),NTIA(2),
     &               NF1AN(2),NF1AA(2),NF1IN(2),NF1II(2),NF1AI(2),
     &               NF2AN(2),NF2AA(2),NF2IN(2),NF2II(2),NF2AI(2)
C
      INTEGER FSDPDAN,FSDPDNA,FSDPDAA,FSDPDIN,FSDPDNI,FSDPDII,FSDPDAI,
     $   FSDPDIA
      COMMON /FSSYMPOP/ FSDPDAN(8,22),FSDPDNA(8,22),FSDPDAA(8,22),
     &                  FSDPDIN(8,22),FSDPDNI(8,22),FSDPDII(8,22),
     &                  FSDPDAI(8,22),FSDPDIA(8,22)
C
      COMMON /ACTIRR/ IRREPOA,IRREPOB,IRREPVA,IRREPVB
C
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA ONE,ONEM,ZILCH/1.D0,-1.D0,0.D0/
C
C     GET T1 ELEMETS FOR TAU
C
C     I0TAIA: T1(AI,A)  ETC..
C
      I0TAIA=MAXCOR-NTAI(1)*IINTFP+1
      I0TAIB=I0TAIA-NTAI(2)*IINTFP
      MXCOR=I0TAIB-1
      CALL FSGETT1(ICORE(I0TAIA),1,90,'AI',9)
      CALL FSGETT1(ICORE(I0TAIB),2,90,'AI',10)
C
      LISTT1=90
      IF(.NOT.LAMBDA) THEN
         LISTAR=63
      ELSE
         LISTAR=146
      ENDIF
      LISTT2=46
C
      NSIZT1A=FSDPDIA(1,9)
      NSIZT1B=FSDPDIA(1,10)
C     
C     I000:  T1A   I010:  T1B
C     
      I000=1
      I010=I000+NSIZT1A*IINTFP
      I020=I010+NSIZT1B*IINTFP
      CALL FSGETT1(ICORE(I000),1,LISTT1,'IA',9)
      CALL FSGETT1(ICORE(I010),2,LISTT1,'IA',10)
C     
      IRREP=DIRPRD(IRREPVA,IRREPVB)
      NUMSYT2=FSDPDII(IRREP,ISYTYP(2,LISTT2))
      DISSYT2=FSDPDAA(IRREP,ISYTYP(1,LISTT2))
      IF(NUMSYT2*DISSYT2.EQ.0) RETURN
C     
C     I000:  T1A  I010: T1B    I020: T2(Nm,iJ)   I030: T2(Nm,Ji)
C     
      I030=I020+NUMSYT2*DISSYT2*IINTFP
      I040=I030+NUMSYT2*DISSYT2*IINTFP
C     
      CALL FSGET(ICORE(I030),1,NUMSYT2,1,IRREP,LISTT2,'AAII')
      CALL FTAU(ICORE(I030),ICORE(I0TAIA),ICORE(I0TAIB),DISSYT2,NUMSYT2,
     $   POPI(1,1),POPI(1,2),VRTA(1,1),VRTA(1,2),IRREP,3,ONE)
      CALL SYMTRA(IRREP,POPI(1,1),POPI(1,2),DISSYT2,ICORE(I030),
     $   ICORE(I020))
C     
C     FIRST TERM
C     
      NUMSYT=FSDPDII(IRREP,ISYTYP(2,LISTAR))
      DISSYT=FSDPDIA(IRREP,ISYTYP(1,LISTAR))
C     
C     I000:  T1A  I010: T1B    I020: T2    I030: T2(TARGET) I040: L2
C     
      I040=I030+NUMSYT*DISSYT*IINTFP
      IF(.NOT.LAMBDA) THEN
         I050=I040
      ELSE
         I050=I040+NUMSYT*DISSYT*IINTFP
      ENDIF
      IF(I050.GT.MAXCOR) CALL INSMEM('OST2IIAI2',I050,MAXCOR)
      IF(NUMSYT*DISSYT.NE.NUMSYT2*DISSYT2*NSIZT1B) THEN
         WRITE(6,*)'@OST2IIAI2-E FATAL ERROR: NUMSYT*DISSYT,',
     $      '.NE.NUMSYT2*DISSYT2*NSIZT1B',NUMSYT,DISSYT,NUMSYT2,
     $      DISSYT2,NSIZT1B
         CALL ERREX
      ENDIF
C     
      IF(.NOT.LAMBDA) THEN
         CALL FSGET(ICORE(I030),1,NUMSYT,1,IRREP,LISTAR,'IAII')
         FACT=ONE
      ELSE
         FACT=ZILCH
      ENDIF
      CALL XGEMM('N','N',NSIZT1B,NUMSYT2*DISSYT2,1,ONEM*W,
     $   ICORE(I010),NSIZT1B,ICORE(I020),1,FACT,ICORE(I030),
     $   NSIZT1B)
      IF(.NOT.LAMBDA) THEN
         CALL FSPUT(ICORE(I030),1,NUMSYT,1,IRREP,LISTAR,'IAII')
      ELSE
         CALL FSGET(ICORE(I040),1,NUMSYT,1,IRREP,LISTAR,'IAII')
         WL=WL+SDOT(NUMSYT*DISSYT,ICORE(I040),1,ICORE(I030),1)
      ENDIF
C     
C     SECOND TERM
C     
      NUMSYT=FSDPDII(IRREP,ISYTYP(2,LISTAR))
      DISSYT=FSDPDAI(IRREP,ISYTYP(1,LISTAR))
C     
C     I000:  T1A  I010: T1B    I020: T2    I030: T2(TARGET)  I040: L2
C     
      I040=I030+NUMSYT*DISSYT*IINTFP
      IF(.NOT.LAMBDA) THEN
         I050=I040
      ELSE
         I050=I040+NUMSYT*DISSYT*IINTFP
      ENDIF
      IF(I050.GT.MAXCOR) CALL INSMEM('OST2IIAI2',I050,MAXCOR)
      IF(NUMSYT*DISSYT.NE.NUMSYT2*DISSYT2*NSIZT1A) THEN
         WRITE(6,*)'@OST2IIAI2-E FATAL ERROR: NUMSYT*DISSYT,',
     $      '.NE.NUMSYT2*DISSYT2*NSIZT1A',NUMSYT,DISSYT,NUMSYT2,
     $      DISSYT2,NSIZT1A
         CALL ERREX
      ENDIF
C     
      IF(.NOT.LAMBDA) THEN
         CALL FSGET(ICORE(I030),1,NUMSYT,1,IRREP,LISTAR,'AIII')
         FACT=ONE
      ELSE
         FACT=ZILCH
      ENDIF
      CALL XGEMM('N','N',NSIZT1A,NUMSYT2*DISSYT2,1,ONEM*W,
     $   ICORE(I000),NSIZT1A,ICORE(I020),1,FACT,ICORE(I030),
     $   NSIZT1A)
      IF(.NOT.LAMBDA) THEN
         CALL FSPUT(ICORE(I030),1,NUMSYT,1,IRREP,LISTAR,'AIII')
      ELSE
         CALL FSGET(ICORE(I040),1,NUMSYT,1,IRREP,LISTAR,'AIII')
         WL=WL+SDOT(NUMSYT*DISSYT,ICORE(I040),1,ICORE(I030),1)
      ENDIF
C     
      RETURN
      END
