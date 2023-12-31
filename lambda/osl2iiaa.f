      SUBROUTINE OSL2IIAA(ICORE,MAXCOR,IUHF,W)
C
C  THIS ROUTINE COMPUTES:
C
C     L2INC(Ab,Mn) = W [ - L2(Ba,Mn) + SUM(i) T1(m,i) * L2(Ba,In) +
C                                    + SUM(I) T1(N,I) * L2(Ba,Mi) +
C                                    - SUM(iJ) TAU(Nm,Ji) * L2(Ba,Ij) ]
C
CEND
C
C CODED PS SEP/93
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYL,DISSYT,POP,VRT,DISTAR
      CHARACTER*4 STRING
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
      COMMON /ACTIRR/ IRREPO(2),IRREPV(2)
C
      COMMON /FLAGS/ IFLAGS(100)
C
      DATA ONE,ONEM,ZILCH,HALF/1.D0,-1.D0,0.D0,0.5D0/
C
      LISTAR=63
      IRREP=DIRPRD(IRREPO(1),IRREPO(2))
      DISTAR=FSDPDII(IRREP,ISYTYP(1,LISTAR))
      NUMTAR=FSDPDAA(IRREP,ISYTYP(2,LISTAR))
      IF(DISTAR*NUMTAR.EQ.0) RETURN
C
C     I000: L2(CONT)
C
         I000=1
         I010=I000+DISTAR*NUMTAR*IINTFP
C
C     FIRST TERM
C
         LISTL=146
         NUMSYL=FSDPDAA(IRREP,ISYTYP(2,LISTL))
         DISSYL=FSDPDII(IRREP,ISYTYP(1,LISTL))
C
         IF(DISTAR.NE.DISSYL) THEN
            WRITE(6,*)'@OSL2IIAA-E FATAL ERROR: DISTAR.NE.DISSYL',
     $           DISTAR,DISSYL
            CALL ERREX
         ENDIF
C
         CALL FSGET(ICORE(I000),1,NUMSYL,1,IRREP,LISTL,'IIAA')
         CALL SSCAL(NUMTAR*DISTAR,ONEM,ICORE(I000),1)
C
C     SECOND AND THIRD TERMS
C
      DO 1000 ISPIN=1,IUHF+1
         NSIZT=NTAI(3-ISPIN)
         LISTL=146
C
         DISSYL=FSDPDII(IRREP,ISYTYP(1,LISTL))
         IF(ISPIN.EQ.1) THEN
            NUMSYL=FSDPDIA(IRREP,ISYTYP(2,LISTL))
            STRING='IIIA'
         ELSE
            NUMSYL=FSDPDAI(IRREP,ISYTYP(2,LISTL))
            STRING='IIAI'
         ENDIF
C
         IF(NSIZT.NE.NUMSYL) THEN
            WRITE(6,*)'@OSL2IIAA-E FATAL ERROR: NSIZT.NE.NUMSYL',
     $           NSIZT,NUMSYL
            CALL ERREX
         ENDIF
C
         IF(DISTAR.NE.DISSYL) THEN
            WRITE(6,*)'@OSL21IIAA-E FATAL ERROR: DISTAR.NE.DISSYL',
     $           DISTAR,DISSYL
            CALL ERREX
         ENDIF
C
C     I000: L2(CONT) I010: T1  I020: L2
C
         I020=I010+NSIZT*IINTFP
         I030=I020+NUMSYL*DISSYL*IINTFP
C
         CALL FSGETT1(ICORE(I010),3-ISPIN,90,'AI',11-ISPIN)
         CALL FSGET(ICORE(I020),1,NUMSYL,1,IRREP,LISTL,STRING)
C
         CALL XGEMM('N','N',DISSYL,1,NUMSYL,ONE,ICORE(I020),DISSYL,
     $        ICORE(I010),NSIZT,ONE,ICORE(I000),DISTAR)
C
1000  CONTINUE
C
C     FOURTH TERM
C
c&test NO CONTRIB IN FIRST TEST
         LISTT=46
         LISTL=46
C
         DISSYL=FSDPDII(IRREP,ISYTYP(1,LISTL))
         NUMSYL=FSDPDII(IRREP,ISYTYP(2,LISTL))
         DISSYT=FSDPDAA(IRREP,ISYTYP(1,LISTT))
         NUMSYT=FSDPDII(IRREP,ISYTYP(2,LISTT))
         NSIZTA=NTAI(1)
         NSIZTB=NTAI(2)
C
C     I000: L2(CONT)  I010: L2   I020:  TAU  I030: T1A   I040: T1B I050: TMP
C
         I020=I010+IINTFP*NUMSYL*DISSYL
         I030=I020+IINTFP*NUMSYT*DISSYT
         I040=I030+IINTFP*NSIZTA
         I050=I040+IINTFP*NSIZTB
         I060=I050+IINTFP*NUMSYT*DISSYT
         IF(I060.GT.MAXCOR) CALL INSMEM('OSL2IIAA',I060,MAXCOR)
C
         CALL FSGETT1(ICORE(I030),1,90,'AI',9)
         CALL FSGETT1(ICORE(I040),2,90,'AI',10)
         CALL FSGET(ICORE(I050),1,NUMSYT,1,IRREP,LISTT,'AAII')
         CALL FTAU(ICORE(I050),ICORE(I030),ICORE(I040),DISSYT,NUMSYT,
     $        POPI(1,1),POPI(1,2),VRTA(1,1),VRTA(1,2),IRREP,3,ONE)
         CALL SYMTRA(IRREP,POPI(1,1),POPI(1,2),DISSYT,
     $        ICORE(I050),ICORE(I020))
         CALL FSGET(ICORE(I010),1,NUMSYL,1,IRREP,LISTL,'IIII')
C
         CALL XGEMM('N','T',DISSYL,1,NUMSYL,ONEM,ICORE(I010),DISSYL,
     $        ICORE(I020),DISSYT,ONE,ICORE(I000),DISTAR)
C
C     SCALE BY W
C
C     I000: L2(TARGET)    I010: L2(CONT)
C
      CALL SYMTRA2(IRREP,VRTI(1,1),VRTI(1,2),DISTAR,NUMTAR,
     $         ICORE(I000),ICORE(I010))
      CALL FSGET(ICORE(I000),1,NUMTAR,1,IRREP,LISTAR,'IIAA')
      CALL SAXPY(NUMTAR*DISTAR,W,ICORE(I010),1,ICORE(I000),1)
      CALL FSPUT(ICORE(I000),1,NUMTAR,1,IRREP,LISTAR,'IIAA')
      RETURN
      END
