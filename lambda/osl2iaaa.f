       SUBROUTINE OSL2IAAA(ICORE,MAXCOR,IUHF,W)
C
C THIS ROUTINE COMPUTES :
C
C     L2INC(Am,Mn) = TERM1A + TERM2A + TERM3A + TERM4A + TERM5A (ISPIN=1)
C
C     L2INC(Na,Mn) = TERM1B + TERM2B + TERM3B + TERM4B + TERM5B (ISPIN=2)
C
C  WITH 
C
C     TERM1A = - L1(a,n)                TERM1B = - L1(A,M)
C     TERM2A = SUM(I) T1(N,I)*L1(a,i)   TERM2B = SUM(i) T1(m,i)*L1(A,I)
C     TERM3A = 1/2 SUM(IJB) T2(BN,IJ)*L2(ba,ij) + SUM(iJb) T2(Nb,Ji)*L2(Ba,Ij)
C     TERM3B = 1/2 SUM(ijb) T2(bm,ij)*L2(BA,IJ) + SUM(IjB) T2(Bm,Ij)*L2(Ab,Ji)
C     TERM4A =   SUM(BI) (T1(B,I)-T1(b,i))*L2(Ba,In) 
C              + SUM(bi) (T1(b,i)-T1(B,I))*L2(ab,ni)
C     TERM4B =   SUM(bi) (T1(b,i)-T1(B,I))*L2(Ab,Mi)
C              + SUM(BI) (T1(B,I)-T1(b,i))*L2(AB,MI)
C     TERM5A = - SUM(Ibj) T1(N,I)*(T1(b,j)-T1(B,J))*L2(ba,ji)
C              - SUM(IBJ) T1(N,I)*(T1(B,J)-T1(b,j))*L2(Ba,Ji)
C     TERM5B = - SUM(iBJ) T1(m,i)*(T1(B,J)-T1(b,j))*L2(BA,JI)
C              - SUM(ijb) T1(m,i)*(T1(b,j)-T1(B,J))*L2(Ab,Ij)
CEND
C
C CODED PS AUG/93
C 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYL,DISSYT,POP,VRT,DISTAR
      CHARACTER*4 STRING,STRTAR
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0(2),DISTAR(2)
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
      DATA ONE,ONEM,ZILCH/1.D0,-1.D0,0.D0/
C
      FACT=ONE
      LISTAR=63
      IRREPT=DIRPRD(IRREPO(1),IRREPO(2))
      NUMTAR=FSDPDAA(IRREPT,ISYTYP(2,LISTAR))
      DISTAR(1)=FSDPDIA(IRREPT,ISYTYP(1,LISTAR))
      DISTAR(2)=FSDPDAI(IRREPT,ISYTYP(1,LISTAR))
      IF((NUMTAR*DISTAR(1).EQ.0).AND.(NUMTAR*DISTAR(2).EQ.0)) RETURN
C
C     READ IN T1 AMPLITUDES REQUIRED.
C     DIRTY TRICKS ARE PLAYED BECAUSE (A,I) IN THE FIRST TERM IS 
C     UNRESTRICTED BUT IN THE SECOND ONE ONLY INACTIVE IS ALLOWED
C
      LISTT1=90
      NSIZTA=FSDPDII(1,9)
      NSIZTB=FSDPDII(1,10)
      IF(NSIZTA.NE.NSIZTB) THEN
         WRITE(6,*)'@OSL2AAIA-E FATAL ERROR: NSIZTA.NE.NSIZTB',
     $      NSIZTA,NSIZTB
         CALL ERREX
      ENDIF
      NSIZTNA=NT(1)
      NSIZTNB=NT(2)
C
C     I0(1):  T1A-T1B  I0(2): T1B-T1A  I010:  T1A  I020:  T1B
C
      I0(1)=1
      I0(2)=I0(1)+MAX(NSIZTA,NSIZTNA)*IINTFP
      I010=I0(2)+MAX(NSIZTB,NSIZTNB)*IINTFP
      I020=I010+NSIZTA*IINTFP
      I030=I020+NSIZTB*IINTFP
      CALL FSGETT1(ICORE(I010),1,LISTT1,'II',9)
      CALL FSGETT1(ICORE(I020),2,LISTT1,'II',10)
      CALL SCOPY(NSIZTA,ICORE(I010),1,ICORE(I0(1)),1)
      CALL SCOPY(NSIZTA,ICORE(I020),1,ICORE(I0(2)),1)
      CALL SAXPY(NSIZTA,ONEM,ICORE(I020),1,ICORE(I0(1)),1)
      CALL SAXPY(NSIZTA,ONEM,ICORE(I010),1,ICORE(I0(2)),1)
      CALL FSPUTT1(ICORE(I0(1)),1,LISTT1,'II',9)
      CALL FSPUTT1(ICORE(I0(2)),2,LISTT1,'II',10)
      CALL GETLST(ICORE(I0(1)),1,1,1,1,LISTT1)
      CALL GETLST(ICORE(I0(2)),1,1,1,2,LISTT1)
      CALL FSPUTT1(ICORE(I010),1,LISTT1,'II',9)
      CALL FSPUTT1(ICORE(I020),2,LISTT1,'II',10)
C
      DO 1000 ISPIN=1,IUHF+1
C
         IF(ISPIN.EQ.1) THEN
            STRTAR='IAAA'
         ELSE
            STRTAR='AIAA'
         ENDIF
         NSIZTAR=NTIA(3-ISPIN)
         IF(NSIZTAR.NE.DISTAR(ISPIN)) THEN
            WRITE(6,*)'@OSL2IAAA-E FATAL ERROR: NSIZTAR.NE.DISTAR',
     $           NSIZTAR,DISTAR(ISPIN)
            CALL ERREX
         ENDIF
         IF(NUMTAR*DISTAR(ISPIN).EQ.0) GOTO 1000
C
C     I0(1):  T1A-T1B  I0(2): T1B-T1A  I000: L2(CONT)
C
         I000=I010
         I010=I000+NSIZTAR*IINTFP
C          
C     FIRST TERM
C
         CALL FSGETT1(ICORE(I000),3-ISPIN,190,'IA',11-ISPIN)
         CALL SSCAL(NSIZTAR,ONEM,ICORE(I000),1)
C
C     SECOND TERM
C
         NSIZL=NTII(3-ISPIN)
         NSIZT=NTAI(ISPIN)
C
C     I000: L2(CONT)  I010:  L1      I020:  T1
C
         I020=I010+NSIZL*IINTFP
         I030=I020+NSIZT*IINTFP
C
         IRREP=IRREPO(3-ISPIN)
         IOFFL=0
         IOFFT=0
         DO 10 IRREPI=1,IRREP-1
            IOFFL=IOFFL+VRTI(IRREPI,3-ISPIN)*POPI(IRREPI,3-ISPIN)*IINTFP
            IOFFT=IOFFT+VRTA(IRREPI,ISPIN)*POPI(IRREPI,ISPIN)*IINTFP
 10      CONTINUE
         CALL FSGETT1(ICORE(I010),3-ISPIN,190,'II',11-ISPIN)
         CALL FSGETT1(ICORE(I020),ISPIN,90,'AI',8+ISPIN)
C
         CALL XGEMM('N','T',VRTI(IRREP,3-ISPIN),1,POPI(IRREP,ISPIN),
     $      ONE,ICORE(I010+IOFFL),VRTI(IRREP,3-ISPIN),ICORE(I020+IOFFT),
     $      1,ONE,ICORE(I000),VRTI(IRREP,3-ISPIN))
C
C     THIRD TERM
C
C     AA AND BB SPIN CASES
C
         LISTT=ISPIN+43
         LISTL=3-ISPIN+143
C     
         DO 100 IRREP=1,NIRREP
C     
            DISSYL=FSDPDII(IRREP,ISYTYP(1,LISTL))
            DISSYT=FSDPDIA(IRREP,ISYTYP(1,LISTT))
            NUMSYL=FSDPDII(IRREP,ISYTYP(2,LISTL))
            NUMSYT=FSDPDII(IRREP,ISYTYP(2,LISTT)) 
C
C     I000: L2(CONT)   I010: L2     I020: T2   I030: TMP
C
            I020=I010+IINTFP*NUMSYL*DISSYL
            I030=I020+IINTFP*NUMSYT*DISSYT
            I040=I030+IINTFP*MAX(NUMSYT*DISSYT,NUMSYL*DISSYL)
            IF(MIN(NUMSYT,NUMSYL,DISSYT,DISSYL).NE.0) THEN
               IF(I040.GT.MAXCOR) CALL INSMEM('OSL2IAAA',I040,MAXCOR)
C  
               CALL OSL2IAA(ICORE(I010),ICORE(I020),ICORE(I030),
     $            ICORE(I000),ISPIN,VRTI(1,3-ISPIN),VRTA(1,ISPIN),
     $            VRTI(1,ISPIN),NSIZTAR,DISSYL,DISSYT,
     &            NUMSYL,NUMSYT,LISTL,LISTT,IRREP,FACT)
            ENDIF      
 100     CONTINUE
C
C       AB SPIN CASE
C
       LISTT=46
       LISTL=146
C
       DO 200 IRREP=1,NIRREP
C
        DISSYL=FSDPDII(IRREP,ISYTYP(1,LISTL))
        NUMSYL=FSDPDII(IRREP,ISYTYP(2,LISTL))
        NUMSYT=FSDPDII(IRREP,ISYTYP(2,LISTT))
        IF(ISPIN.EQ.1) THEN
           DISSYT=FSDPDAI(IRREP,ISYTYP(1,LISTT))
           STRING='AIII'
        ELSE
           DISSYT=FSDPDIA(IRREP,ISYTYP(1,LISTT))
           STRING='IAII'
        ENDIF
C
C     I000: L2(CONT)  I010: L2   I020: T2   I030: TMP   I040: SCR
C
        I020=I010+IINTFP*NUMSYL*DISSYL
        I030=I020+IINTFP*NUMSYT*DISSYT
        I040=I030+IINTFP*NUMSYT*DISSYT
        IF(MIN(NUMSYT,NUMSYL,DISSYT,DISSYL).NE.0) THEN
           I050=I040+3*IINTFP*MAX(DISSYT,DISSYL,NUMSYL,NUMSYT)
           IF(I050.GT.MAXCOR) CALL INSMEM('OSL2IAAA',I050,MAXCOR)
C     
           CALL OSL2IAB(ICORE(I010),ICORE(I020),ICORE(I030),
     $        ICORE(I000),ISPIN,VRTA(1,ISPIN),VRTI(1,3-ISPIN),
     $        VRTI(1,ISPIN),POPI(1,1),POPI(1,2),NSIZTAR,DISSYL,DISSYT,
     $        NUMSYL,NUMSYT,LISTL,LISTT,IRREP,ICORE(I040),FACT,STRING)
        ENDIF
 200  CONTINUE
C
C     FOURTH TERM
C     WE ONLY NEED TO CONSIDER ONE DPD IRREP -- THE TOTALLY SYMMETRIC ONE.
C
      LISTL1=136-ISPIN
      LISTL2=138-ISPIN
      NLDSZ1=IRPDPD(1,ISYTYP(1,LISTL1))
      NLDSZ2=IRPDPD(1,ISYTYP(1,LISTL2))
      NLDIS1=FSDPDIA(1,ISYTYP(2,LISTL1))
      NLDIS2=FSDPDIA(1,ISYTYP(2,LISTL2))
C
C     I0(1):  T1A-T1B  I0(2): T1B-T1A  I000: L2(CONT) I010: L2
C
      I020=I010+NLDSZ2*NLDIS2*IINTFP
      IF(I020.GT.MAXCOR)CALL INSMEM('OSL2IAAA',I020,MAXCOR)
C
      CALL FSGET(ICORE(I010),1,NLDIS2,1,1,LISTL2,'NNIA')
      CALL XGEMM('N','N',1,NSIZTAR,NLDSZ2,ONE,ICORE(I0(ISPIN)),
     &           1,ICORE(I010),NLDSZ2,ONE,ICORE(I000),1)
C
      I020=I010+NLDSZ1*NLDIS1*IINTFP
      IF(I020.GT.MAXCOR)CALL INSMEM('OSL2IAAA',I020,MAXCOR)
C
      CALL FSGET(ICORE(I010),1,NLDIS1,1,1,LISTL1,'NNIA')
      CALL XGEMM('N','N',1,NSIZTAR,NLDSZ1,ONEM,ICORE(I0(3-ISPIN)),
     &           1,ICORE(I010),NLDSZ1,ONE,ICORE(I000),1)
C
C     FIFTH TERM
C
C     AA AND BB SPIN CASES
C
         LISTT=46
         LISTL=3-ISPIN+143
C     
         DO 300 IRREP=1,NIRREP
C     
            DISSYL=FSDPDNI(IRREP,ISYTYP(1,LISTL))
            NUMSYL=FSDPDNI(IRREP,ISYTYP(2,LISTL))
            IF(ISPIN.EQ.1) THEN
               NUMSYT=FSDPDIN(IRREP,ISYTYP(2,LISTT)) 
               DISSYT=FSDPDAN(IRREP,ISYTYP(1,LISTT))
            ELSE
               NUMSYT=FSDPDNI(IRREP,ISYTYP(2,LISTT)) 
               DISSYT=FSDPDNA(IRREP,ISYTYP(1,LISTT))
            ENDIF
            NSIZT=NTAI(ISPIN)
C
C     I000: L2(CONT)   I010: L2     I020: T2   I030: TMP   I040: T1
C
            I020=I010+IINTFP*NUMSYL*DISSYL
            I030=I020+IINTFP*NUMSYT*DISSYT
            I040=I030+IINTFP*MAX(NUMSYT*DISSYT,NUMSYL*DISSYL)
            I050=I040+IINTFP*NSIZT
            IF(MIN(NUMSYT,NUMSYL,DISSYT,DISSYL).NE.0) THEN
               IF(I050.GT.MAXCOR) CALL INSMEM('OSL2IAAA',I050,MAXCOR)
C
               CALL FSGETT1(ICORE(I040),ISPIN,90,'AI',8+ISPIN)
C
               CALL OSL2IAC(ICORE(I010),ICORE(I020),ICORE(I030),
     $            ICORE(I000),ICORE(I0(3-ISPIN)),ICORE(I040),
     $            ISPIN,VRTI(1,3-ISPIN),VRT(1,3-ISPIN),
     $            VRTA(1,ISPIN),POP(1,3-ISPIN),POPI(1,ISPIN),
     $            NSIZTAR,DISSYL,DISSYT,
     $            NUMSYL,NUMSYT,LISTL,IRREP,FACT)
            ENDIF      
 300     CONTINUE
C
C       AB SPIN CASE
C
       LISTT=43+ISPIN
       LISTL=146
C
       DO 400 IRREP=1,NIRREP
C
        IF(ISPIN.EQ.1) THEN
           DISSYL=FSDPDNI(IRREP,ISYTYP(1,LISTL))
           NUMSYL=FSDPDNI(IRREP,ISYTYP(2,LISTL))
           STRING='NINI'
        ELSE
           DISSYL=FSDPDIN(IRREP,ISYTYP(1,LISTL))
           NUMSYL=FSDPDIN(IRREP,ISYTYP(2,LISTL))
           STRING='ININ'
        ENDIF
        NUMSYT=FSDPDNI(IRREP,ISYTYP(2,LISTT))
        DISSYT=FSDPDNA(IRREP,ISYTYP(1,LISTT))
        NSIZT=NTAI(ISPIN)
C
C     I000: L2(CONT)  I010: L2 I020: T2 I030: T1 I040: TMP 
C
        I020=I010+IINTFP*NUMSYL*DISSYL
        I030=I020+IINTFP*NUMSYT*DISSYT
        I040=I030+IINTFP*NSIZT
        IF(MIN(NUMSYT,NUMSYL,DISSYT,DISSYL).NE.0) THEN
           I050=I040+IINTFP*MAX(DISSYT*NUMSYT,DISSYL*NUMSYL)
           IF(I050.GT.MAXCOR) CALL INSMEM('OSL2IAAA',I050,MAXCOR)
C
           CALL FSGETT1(ICORE(I030),ISPIN,90,'AI',8+ISPIN)
C
           CALL OSL2IAD(ICORE(I010),ICORE(I020),ICORE(I040),
     $        ICORE(I000),ICORE(I0(ISPIN)),ICORE(I030),
     $        ISPIN,VRT(1,ISPIN),VRTA(1,ISPIN),
     $        VRTI(1,3-ISPIN),VRT(1,ISPIN),POP(1,ISPIN),POPI(1,ISPIN),
     $        NSIZTAR,DISSYL,DISSYT,
     $        NUMSYL,NUMSYT,LISTL,LISTT,IRREP,ICORE(I040),FACT,STRING)
        ENDIF
 400  CONTINUE
C
C     SCALE BY W
C
C     I000: L2(CONT)    I010: L2(TARGET)
C
      CALL FSGET(ICORE(I010),1,NUMTAR,1,IRREPT,LISTAR,STRTAR)
      CALL SAXPY(NUMTAR*DISTAR(ISPIN),W,ICORE(I000),1,ICORE(I010),1)
      CALL FSPUT(ICORE(I010),1,NUMTAR,1,IRREPT,LISTAR,STRTAR)
C
1000  CONTINUE
      RETURN
      END
