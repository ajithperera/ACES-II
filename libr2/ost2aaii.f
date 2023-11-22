      SUBROUTINE OST2AAII(ICORE,MAXCOR,IUHF,W,LAMBDA,WL)
C
C
C     THIS ROUTINE CALCULATE ACTIVE,ACTIVE-INACTIVE,INACTIVE T2 CONTRIBUTION
C
C     T2(Mn,Ab)= -W * (T2(nM,aB)+T1(a,n)*T1(B,M)) = -W * TAU(nM,aB)   
C
CEND
C
CPROGRAMED PS NOV/92
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,POP,VRT
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
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA ONE,ONEM/1.D0,-1.D0/
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
      DO 100 IRREP=1,NIRREP
C     
         NUMSYT=FSDPDAA(IRREP,ISYTYP(2,LISTT2))
         IF(NUMSYT.EQ.0) GOTO 100
C     
         DISSYT=FSDPDII(IRREP,ISYTYP(1,LISTT2))
         STRING='IIAA'
         IF(DISSYT.EQ.0) GOTO 100
C     
         IF(NUMSYT.NE.1) THEN
            WRITE(6,*)'@OST2IIAA-E FATAL ERROR: NUMSYT,',
     $         NUMSYT
            CALL ERREX
         ENDIF
C     
C     I000:  T1A  I010:  T1B  I020: T2(TARGET)   I030: T2
C     
         I030=I020+NUMSYT*DISSYT*IINTFP
         I040=I030+NUMSYT*DISSYT*IINTFP
C     
         CALL FSGET(ICORE(I020),1,NUMSYT,1,IRREP,LISTT2,STRING)
         CALL SYMTRA2(IRREP,VRTI(1,1),VRTI(1,2),DISSYT,
     $      NUMSYT,ICORE(I020),ICORE(I030))
         CALL FTAU(ICORE(I030),ICORE(I010),ICORE(I000),DISSYT,
     $      NUMSYT,POPA(1,2),POPA(1,1),VRTI(1,2),VRTI(1,1),
     $      IRREP,3,ONE)
C
         CALL FSGET(ICORE(I020),1,NUMSYT,1,IRREP,LISTAR,STRING)
         IF(.NOT.LAMBDA) THEN
            CALL SAXPY(NUMSYT*DISSYT,ONEM*W,ICORE(I030),1,ICORE(I020),1)
            CALL FSPUT(ICORE(I020),1,NUMSYT,1,IRREP,LISTAR,STRING)
         ELSE
            WL=WL+ONEM*W*SDOT(NUMSYT*DISSYT,ICORE(I020),1,ICORE(I030),1)
         ENDIF
C     
 100  CONTINUE
      RETURN
      END