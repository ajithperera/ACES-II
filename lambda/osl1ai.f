      

      SUBROUTINE OSL1AI(ICORE,MAXCOR,IUHF,W)
C
C  THIS ROUTINE COMPUTES:
C
C     L1INC(N,I)= W [ SUM(X) TERMXA ]  (ISPIN=1)
C
C     L1INC(m,i)= W [ SUM(X) TERMXB ]  (ISPIN=2)
C
C  WITH
C
C     TERM1A  = L2(Nm,Mi)                    TERM1B = L2(Nm,In)
C     TERM2A  = + SUM(A) T1(A,M) * L2(Na,Mi) + SUM(a) T1(a,n) * L2(Am,Mi)
C     TERM2B  = + SUM(a) T1(a,n) * L2(Am,In) + SUM(A) T1(A,M) * L2(Na,In)
C     TERM3A  = - SUM(j) T1(m,j) * L2(Nm,Ji)
C     TERM3B  = - SUM(J) T1(N,J) * L2(Nm,Ij)
C     TERM4A  = + SUM(A) T2(Am,Mn) * L1(a,i)
C     TERM4B  = + SUM(a) T2(Na,Mn) * L1(A,I)
C     TERM5A  = + SUM(Ab) TAU(Ba,Mn) * L2(Ab,Mi)
C     TERM5B  = + SUM(aB) TAU(Ab,Mn) * L2(Ba,In)
C     TERM6A  = - SUM(AJ)  T2(Am,Jn) * L2(am,ij)
C     TERM6B  = - SUM(aj)  T2(Na,Mj) * L2(AN,IJ)
C     TERM7A  = - SUM(aj) TAU(am,nj) * L2(Am,Ji)
C               - SUM(jA) TAU(Am,Mj) * L2(Na,Ji)
C     TERM7B  = - SUM(AJ) TAU(AN,MJ) * L2(Na,Ij)
C               - SUM(Ja) TAU(Na,Jn) * L2(Am,Ij)
C     TERM8A  = - SUM(jAb) T1(m,j) * TAU(Ab,Mn) * L2(Ba,Ji)
C     TERM8B  = - SUM(JaB) T1(N,J) * TAU(Ba,Mn) * L2(Ab,Ij)
C     TERM9A  = - SUM(jAb) T2(Am,Mn) * (T1(b,j)-T1(B,J)) * L2(ab,ij)
C               - SUM(JAB) T2(Am,Mn) * (T1(B,J)-T1(b,j)) * L2(Ba,Ji)
C     TERM9B  = - SUM(JaB) T2(Na,Mn) * (T1(B,J)-T1(b,j)) * L2(AB,IJ)
C               - SUM(jab) T2(Na,Mn) * (T1(b,j)-T1(B,J)) * L2(Ab,Ij)
C     TERM10A = + SUM(JAB) T1(A,M) * T2(Bm,Jn) * L2(ab,ij)
C               + SUM(jAb) T1(A,M) * T2(bm,jn) * L2(Ba,Ji)
C               - SUM(jaB) T1(a,n) * T2(Bm,Mj) * L2(Ab,Ji)
C     TERM10B = + SUM(jab) T1(a,n) * T2(Nb,Mj) * L2(AB,IJ)
C               + SUM(JAb) T1(a,n) * T2(BN,JM) * L2(Ab,Ij)
C               - SUM(JAb) T1(A,M) * T2(Nb,Jn) * L2(Ba,Ij)
C
CEND
C
C CODED PS AUG/93
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYL,DISSYT,POP,VRT,DISSYT2,DISSYL2
      CHARACTER*4 STRING
C
      DIMENSION ICORE(MAXCOR)
      DIMENSION I0(2),NSTAR(2)
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
      FACT=ONE
      NSTAR(1)=NTAI(1)
      NSTAR(2)=NTAI(2)
      IF((NSTAR(1).EQ.0).AND.(NSTAR(2).EQ.0)) RETURN
C
C     READ IN T1 AMPLITUDES REQUIRED.
C     DIRTY TRICKS ARE PLAYED BECAUSE (A,I) IN THE FIRST TERM IS 
C     UNRESTRICTED BUT IN THE SECOND ONE ONLY INACTIVE IS ALLOWED
C
      LISTT1=90
      NSIZTA=FSDPDII(1,9)
      NSIZTB=FSDPDII(1,10)
      IF(NSIZTA.NE.NSIZTB) THEN
         WRITE(6,*)'@OSL1AI-E FATAL ERROR: NSIZTA.NE.NSIZTB',
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
      I000=I010
      DO 1000 ISPIN=1,IUHF+1
C
         NSIZTAR=NSTAR(ISPIN)
         IF(NSIZTAR.EQ.0) GOTO 1000
C
C     I0(1):  T1A-T1B  I0(2): T1B-T1A  I000: L1(CONT)
C
         I010=I000+NSIZTAR*IINTFP
C          
C     FIRST TERM
C
         LISTL=146
         IRREP=DIRPRD(IRREPV(1),IRREPV(2))
         DISSYL=FSDPDAA(IRREP,ISYTYP(1,LISTL))
         IF(ISPIN.EQ.1) THEN
            NUMSYL=FSDPDAI(IRREP,ISYTYP(2,LISTL))
            STRING='AAAI'
         ELSE
            NUMSYL=FSDPDIA(IRREP,ISYTYP(2,LISTL))
            STRING='AAIA'
         ENDIF
         IF(NSIZTAR.NE.NUMSYL) THEN
            WRITE(6,*)'@OSL21AI-E FATAL ERROR: NSIZTAR.NE.NUMSYL',
     $           NSIZTAR,NUMSYL
            CALL ERREX
         ENDIF
         CALL FSGET(ICORE(I000),1,NUMSYL,1,IRREP,LISTL,STRING)
C
C     SECOND AND TENTH TERM
C     (WITH L2'(Na,Mi)=L2(Na,Mi)+SUM(JB) T2(Bm,Jn)*L2(ab,ij)
C                               +SUM(jb) T2(bm,jn)*L2(Ba,Ji)
C           L2'(Am,Mi)=L2(Am,Mi)-SUM(jB) T2(Bm,Mj)*L2(Ba,Ji)      ISPIN=1
C           L2'(Am,In)=L2(Am,In)+SUM(jb) T2(Nb,Mj)*L2(AB,IJ)
C                               +SUM(jB) T2(BN,JM)*L2(Ab,Ij)
C           L2'(Na,In)=L2(Na,In)-SUM(Jb) T2(Nb,Jn)*L2(Ba,Ij)      ISPIN=2
C
C     AA AND BB SPINCASE
C
         NSIZT=NTIA(ISPIN)
         LISTL=146
         IRREP=DIRPRD(IRREPO(ISPIN),IRREPV(ISPIN))
C
         IF(ISPIN.EQ.1) THEN
            DISSYL=FSDPDAI(IRREP,ISYTYP(1,LISTL))
            NUMSYL=FSDPDAI(IRREP,ISYTYP(2,LISTL))
            STRING='AIAI'
         ELSE
            DISSYL=FSDPDIA(IRREP,ISYTYP(1,LISTL))
            NUMSYL=FSDPDIA(IRREP,ISYTYP(2,LISTL))
            STRING='IAIA'
         ENDIF
C
         IF(NSIZT.NE.DISSYL) THEN
            WRITE(6,*)'@OSL21AI-E FATAL ERROR: NSIZT.NE.DISSYL',
     $           NSIZT,DISSYL
            CALL ERREX
         ENDIF
C
         IF(NSIZTAR.NE.NUMSYL) THEN
            WRITE(6,*)'@OSL21AI-E FATAL ERROR: NSIZTAR.NE.NUMSYL',
     $           NSIZTAR,NUMSYL
            CALL ERREX
         ENDIF
C
C     I000: L1(CONT) I010: T1  I020: L1
C
         I020=I010+NSIZT*IINTFP
         I030=I020+NUMSYL*DISSYL*IINTFP
C
         CALL FSGETT1(ICORE(I010),ISPIN,90,'IA',8+ISPIN)
         CALL FSGET(ICORE(I020),1,NUMSYL,1,IRREP,LISTL,STRING)
C
C           NOW THE CORRECTION FOR TENTH TERM
C
C
C           AA AND BB CASES
C
            LISTT2=35+ISPIN
            LISTL2=136-ISPIN
            NUMSYL2=FSDPDII(IRREP,ISYTYP(2,LISTL2))
            DISSYL2=FSDPDII(IRREP,ISYTYP(1,LISTL2))
            NUMSYT2=FSDPDII(IRREP,ISYTYP(2,LISTT2))
            DISSYT2=FSDPDAA(IRREP,ISYTYP(1,LISTT2))
C
C     I000: L1(CONT) I010: L1  I020: L2  I030: T2(2)   I040:  L2(2)
C
            I040=I030+NUMSYT2*DISSYT2*IINTFP
            I050=I040+NUMSYL2*DISSYL2*IINTFP
C
            CALL FSGET(ICORE(I040),1,NUMSYL2,1,IRREP,LISTL2,'IIII')
            CALL FSGET(ICORE(I030),1,NUMSYT2,1,IRREP,LISTT2,'AAII')
C
            CALL XGEMM('N','N',1,NUMSYL2,DISSYL2,ONEM,ICORE(I030),
     $        1,ICORE(I040),DISSYL2,ONE,ICORE(I020),1)
C
C           AB CASE
C
            LISTT2=36-ISPIN
            LISTL2=138-ISPIN
            NUMSYL2=FSDPDII(IRREP,ISYTYP(2,LISTL2))
            DISSYL2=FSDPDII(IRREP,ISYTYP(1,LISTL2))
            NUMSYT2=FSDPDII(IRREP,ISYTYP(2,LISTT2))
            DISSYT2=FSDPDAA(IRREP,ISYTYP(1,LISTT2))
C
C     I000: L1(CONT) I010: L1  I020: L2  I030: T2(2)   I040:  L2(2)
C
            I040=I030+NUMSYT2*DISSYT2*IINTFP
            I050=I040+NUMSYL2*DISSYL2*IINTFP
C
            CALL FSGET(ICORE(I040),1,NUMSYL2,1,IRREP,LISTL2,'IIII')
            CALL FSGET(ICORE(I030),1,NUMSYT2,1,IRREP,LISTT2,'AAII')
C
            CALL XGEMM('N','N',1,NUMSYL2,DISSYL2,ONEM,ICORE(I030),
     $        1,ICORE(I040),DISSYL2,ONE,ICORE(I020),1)
C
C           END OF TENTH CORRECTION
C
C
C     I000: L1(CONT) I010: L1  I020: L2 
C
         CALL XGEMM('T','N',1,NUMSYL,DISSYL,ONE,ICORE(I010),NSIZT,
     $        ICORE(I020),DISSYL,ONE,ICORE(I000),1)
C
C     AB SPINCASE
C
         NSIZT=NTIA(3-ISPIN)
         LISTL=146
         IRREP=DIRPRD(IRREPO(ISPIN),IRREPV(ISPIN))
C
         IF(ISPIN.EQ.1) THEN
            DISSYL=FSDPDIA(IRREP,ISYTYP(1,LISTL))
            NUMSYL=FSDPDAI(IRREP,ISYTYP(2,LISTL))
            STRING='IAAI'
         ELSE
            DISSYL=FSDPDAI(IRREP,ISYTYP(1,LISTL))
            NUMSYL=FSDPDIA(IRREP,ISYTYP(2,LISTL))
            STRING='AIIA'
         ENDIF
C
         IF(NSIZT.NE.DISSYL) THEN
            WRITE(6,*)'@OSL21AI-E FATAL ERROR: NSIZT.NE.DISSYL',
     $           NSIZT,DISSYL
            CALL ERREX
         ENDIF
C
         IF(NSIZTAR.NE.NUMSYL) THEN
            WRITE(6,*)'@OSL21AI-E FATAL ERROR: NSIZTAR.NE.NUMSYL',
     $           NSIZTAR,NUMSYL
            CALL ERREX
         ENDIF
C
C     I000: L1(CONT) I010: T1  I020: L1
C
         I020=I010+NSIZT*IINTFP
         I030=I020+NUMSYL*DISSYL*IINTFP
C
         CALL FSGETT1(ICORE(I010),3-ISPIN,90,'IA',11-ISPIN)
         CALL FSGET(ICORE(I020),1,NUMSYL,1,IRREP,LISTL,STRING)
C
C           NOW THE CORRECTION FOR TENTH TERM
C
            LISTT2=37+ISPIN
            LISTL2=137+ISPIN
            NUMSYL2=FSDPDII(1,ISYTYP(2,LISTL2))
            DISSYL2=FSDPDII(1,ISYTYP(1,LISTL2))
            NUMSYT2=FSDPDII(1,ISYTYP(2,LISTT2))
            DISSYT2=FSDPDAA(1,ISYTYP(1,LISTT2))
C
C     I000: L1(CONT) I010: L1  I020: L2  I030: T2(2)   I040:  L2(2)
C
            I040=I030+NUMSYT2*DISSYT2*IINTFP
            I050=I040+NUMSYL2*DISSYL2*IINTFP
C
            CALL FSGET(ICORE(I040),1,NUMSYL2,1,1,LISTL2,'IIII')
            CALL FSGET(ICORE(I030),1,NUMSYT2,1,1,LISTT2,'AAII')
C
            IOFF=I040
            IRREPI=IRREPV(ISPIN)
            DO 200 IRREP1=1,IRREPI-1
              IOFF=IOFF+POPI(IRREP1,3-ISPIN)*VRTI(IRREP1,3-ISPIN)*
     $              DISSYL2*IINTFP
 200        CONTINUE
            CALL XGEMM('N','N',1,NUMSYL*DISSYL,DISSYL2,ONEM,ICORE(I030),
     $        1,ICORE(IOFF),DISSYL2,ONE,ICORE(I020),1)
C
C           END OF TENTH CORRECTION
C
C
C     I000: L1(CONT) I010: L1  I020: L2
C
         CALL XGEMM('T','N',1,NUMSYL,DISSYL,ONE,ICORE(I010),NSIZT,
     $        ICORE(I020),DISSYL,ONE,ICORE(I000),1)
C
C     THIRD AND EIGHTH TERM
C     (WITH L2'(Nm,Ji)=L2(Nm,Ji)+ SUM(Ab) TAU(Ab,Mn)*L2(Ba,Ji)
C      AND  L2'(Nm,Ij)=l2(Nm,Ij)+ SUM(aB) TAU(Ba,Mn)*L2(Ab,Ij)  )
C
         NSIZT=NTAI(3-ISPIN)
         LISTL=146
         IRREP=DIRPRD(IRREPV(3-ISPIN),IRREPV(ISPIN))
C
         DISSYL=FSDPDAA(IRREP,ISYTYP(1,LISTL))
         NUMSYL=FSDPDII(IRREP,ISYTYP(2,LISTL))
         STRING='AAII'
C
C     I000: L1(CONT) I010: T1  I020: L2
C
         I020=I010+NSIZT*IINTFP
         I030=I020+NUMSYL*DISSYL*IINTFP
C
         CALL FSGETT1(ICORE(I010),3-ISPIN,90,'AI',11-ISPIN)
         CALL FSGET(ICORE(I020),1,NUMSYL,1,IRREP,LISTL,STRING)
C
C          NOW CORRECT FOR THE EIGHTH TERM
C
           LISTT2=46
           LISTL2=146
           DISSYT2=FSDPDII(IRREP,ISYTYP(1,LISTT2))
           NUMSYT2=FSDPDAA(IRREP,ISYTYP(2,LISTT2))
           DISSYL2=FSDPDII(IRREP,ISYTYP(1,LISTL2))
           NUMSYL2=FSDPDII(IRREP,ISYTYP(2,LISTL2))
           NSIZTA=NTIA(1)
           NSIZTB=NTIA(2)
C
C  I000: L1(CONT) I010: T1  I020: L2 I030: TAU  I040: T1A I050: T1B I060: TMP
C
           I040=I030+NUMSYT2*DISSYT2*IINTFP
           I050=I040+NSIZTA*IINTFP
           I060=I050+NSIZTB*IINTFP
           I070=I060+NUMSYT2*DISSYT2*IINTFP
           IF(I070.GT.MAXCOR) CALL INSMEM('OSL1AI',I070,MAXCOR)
C
           CALL FSGETT1(ICORE(I040),1,90,'IA',9)
           CALL FSGETT1(ICORE(I050),2,90,'IA',10)
           CALL FSGET(ICORE(I060),1,NUMSYT2,1,IRREP,LISTT2,'IIAA')
           CALL FTAU(ICORE(I060),ICORE(I040),ICORE(I050),DISSYT2,NUMSYT2
     $        ,POPA(1,1),POPA(1,2),VRTI(1,1),VRTI(1,2),IRREP,3,ONE)
           CALL SYMTRA2(IRREP,VRTI(1,1),VRTI(1,2),DISSYT2,NUMSYT2,
     $        ICORE(I060),ICORE(I030))
C
C  I000: L1(CONT) I010: T1  I020: L2 I030: TAU  I040: L2(2)
C
           I050=I040+NUMSYL2*DISSYL2*IINTFP
           IF(I050.GT.MAXCOR) CALL INSMEM('OSL1AI',I050,MAXCOR)
C
           CALL FSGET(ICORE(I040),1,NUMSYL2,1,IRREP,LISTL2,'IIII')
C
           CALL XGEMM('T','N',1,NUMSYL2,DISSYL2,ONE,ICORE(I030),DISSYT2,
     $        ICORE(I040),DISSYL2,ONE,ICORE(I020),1)
C
C          END OF EIGTH TERM CORRECTION
C
C
C     I000: L1(CONT) I010: T1  I020: L2  I030,I040,I050: SCR
C
         IF(ISPIN.EQ.1) THEN
            I040=I030+3*DISSYL*IINTFP
            I050=I040+3*DISSYL*IINTFP
            I060=I050+3*DISSYL*IINTFP
            CALL SYMTR1(IRREP,POPI(1,ISPIN),POPI(1,3-ISPIN),DISSYL,
     $           ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050))
         ENDIF
         IRREPJ=IRREPV(3-ISPIN)
         IOFF=I020
         DO 300 IRREPL=1,IRREPJ-1
            IRREPI=DIRPRD(IRREP,IRREPL)
            IOFF=IOFF+POPI(IRREPL,ISPIN)*POPI(IRREPI,3-ISPIN)*IINTFP
 300     CONTINUE
C
         CALL XGEMM('N','T',1,NSIZTAR,NSIZT,ONEM,ICORE(I010),1,
     $        ICORE(IOFF),NSIZTAR,ONE,ICORE(I000),1)
C
C     FOURTH AND NINTH TERMS
C     (WITH: L1'(a,i)=L1(a,i)+SUM(bj) (T1(b,j)-T1(B,J))*L2(ba,ij)
C                            -SUM(Bj) (T1(B,J)-T1(b,j))*L2(Ba,Ji)
C      AND:  L1'(A,I)=L1(A,I)+SUM(JB) (T1(B,J)-T1(b,j))*L2(BA,IJ)
C                            -SUM(jb) (T1(b,j)-T1(B,J))*L2(Ab,Ij)
C      NOTE:  B,J,b,j ARE OF TYPE 'N'  )
C
         NSIZL=NTII(3-ISPIN)
         LISTT=46
         IRREP=DIRPRD(IRREPO(ISPIN),IRREPO(3-ISPIN))
C
         NUMSYT=FSDPDAA(IRREP,ISYTYP(2,LISTT))
         IF(ISPIN.EQ.1) THEN
            DISSYT=FSDPDIA(IRREP,ISYTYP(1,LISTT))
            STRING='IAAA'
         ELSE
            DISSYT=FSDPDAI(IRREP,ISYTYP(1,LISTT))
            STRING='AIAA'
         ENDIF
C
C     I000: L1(CONT) I010: L1  
C
         I020=I010+NSIZL*IINTFP
C
         CALL FSGETT1(ICORE(I010),3-ISPIN,190,'II',11-ISPIN)
C
C           NOW CORRECTION FOR NINTH TERM 
C                 (ONLY TOTALLY SYMMETRIC IRREP REQUIRED)
C
C  AA AND BB SPINCASES
C
            LISTL2=136-ISPIN
            NUMSYL2=FSDPDII(1,ISYTYP(2,LISTL2))
            DISSYL2=IRPDPD(1,ISYTYP(1,LISTL2))
C
C     I0(1):  T1A-T1B  I0(2): T1B-T1A  I000: L1(CONT) I010: L1  I020: L2   
C
            I030=I020+NUMSYL2*DISSYL2*IINTFP
C
            CALL FSGET(ICORE(I020),1,NUMSYL2,1,1,LISTL2,'NNII')
C
            CALL XGEMM('N','N',1,NUMSYL2,DISSYL2,ONE,ICORE(I0(3-ISPIN))
     $        ,1,ICORE(I020),DISSYL2,ONE,ICORE(I010),1)
C
C  AB SPINCASE
C
            LISTL2=138-ISPIN
            NUMSYL2=FSDPDII(1,ISYTYP(2,LISTL2))
            DISSYL2=IRPDPD(1,ISYTYP(1,LISTL2))
C
C     I0(1):  T1A-T1B  I0(2): T1B-T1A  I000: L1(CONT) I010: L1  I020: L2   
C
            I030=I020+NUMSYL2*DISSYL2*IINTFP
C
            CALL FSGET(ICORE(I020),1,NUMSYL2,1,1,LISTL2,'NNII')
C
            CALL XGEMM('N','N',1,NUMSYL2,DISSYL2,ONEM,ICORE(I0(ISPIN)),
     $        1,ICORE(I020),DISSYL2,ONE,ICORE(I010),1)
C
C           END OF NINTH TERM CONTRIBUTION
C
C
C     I000: L1(CONT) I010: L1  I020: T2
C
         I030=I020+NUMSYT*DISSYT*IINTFP
C
         CALL FSGET(ICORE(I020),1,NUMSYT,1,IRREP,LISTT,STRING)
C
         IRREPC=IRREPV(ISPIN)
         IOFF=I010
         DO 400 IRREPL=1,IRREPC-1
            IOFF=IOFF+POPI(IRREPL,3-ISPIN)*VRTI(IRREPL,3-ISPIN)*IINTFP
 400     CONTINUE
C
         CALL XGEMM('T','N',1,NSIZTAR,DISSYT,ONE,ICORE(I020),DISSYT,
     $        ICORE(IOFF),DISSYT,ONE,ICORE(I000),1)
C
C     FIFTH TERM
C
         LISTT=46
         LISTL=146
C
         IRREP=DIRPRD(IRREPO(ISPIN),IRREPO(3-ISPIN))
C
         DISSYL=FSDPDII(IRREP,ISYTYP(1,LISTL))
         DISSYT=FSDPDII(IRREP,ISYTYP(1,LISTT))
         NUMSYT=FSDPDAA(IRREP,ISYTYP(2,LISTT))
         IF(ISPIN.EQ.1) THEN
            NUMSYL=FSDPDAI(IRREP,ISYTYP(2,LISTL))
            STRING='IIAI'
         ELSE
            NUMSYL=FSDPDIA(IRREP,ISYTYP(2,LISTL))
            STRING='IIIA'
         ENDIF
         NSIZTA=NTIA(1)
         NSIZTB=NTIA(2)
C
C     I000: L2(CONT)  I010: L2   I020:  TAU  I030: T1A   I040: T1B I050: TMP
C
         I020=I010+IINTFP*NUMSYL*DISSYL
         I030=I020+IINTFP*NUMSYT*DISSYT
         I040=I030+IINTFP*NSIZTA
         I050=I040+IINTFP*NSIZTB
         I060=I050+IINTFP*NUMSYT*DISSYT
         IF(I060.GT.MAXCOR) CALL INSMEM('OSL1AI',I060,MAXCOR)
C
         CALL FSGETT1(ICORE(I030),1,90,'IA',9)
         CALL FSGETT1(ICORE(I040),2,90,'IA',10)
         CALL FSGET(ICORE(I050),1,NUMSYT,1,IRREP,LISTT,'IIAA')
         CALL FTAU(ICORE(I050),ICORE(I030),ICORE(I040),DISSYT,NUMSYT,
     $        POPA(1,1),POPA(1,2),VRTI(1,1),VRTI(1,2),IRREP,3,ONE)
         CALL SYMTRA2(IRREP,VRTI(1,1),VRTI(1,2),DISSYT,NUMSYT,
     $        ICORE(I050),ICORE(I020))
         CALL FSGET(ICORE(I010),1,NUMSYL,1,IRREP,LISTL,STRING)
C
         CALL XGEMM('T','N',1,NUMSYL,DISSYL,ONE,ICORE(I020),DISSYT,
     $        ICORE(I010),DISSYL,ONE,ICORE(I000),1)
C
C     SIXTH TERM
C
         LISTT=35+ISPIN
         LISTL=136-ISPIN
C
         IRREP=DIRPRD(IRREPO(3-ISPIN),IRREPV(3-ISPIN))
C
         DISSYL=FSDPDII(IRREP,ISYTYP(1,LISTL))
         NUMSYT=FSDPDII(IRREP,ISYTYP(2,LISTT))
         DISSYT=FSDPDAA(IRREP,ISYTYP(1,LISTT))
         NUMSYL=FSDPDAI(IRREP,ISYTYP(2,LISTL))
C
C     I000: L2(CONT)  I010: L2   I020:  T2
C
         I020=I010+IINTFP*NUMSYL*DISSYL
         I030=I020+IINTFP*NUMSYT*DISSYT
         IF(I030.GT.MAXCOR) CALL INSMEM('OSL1AI',I030,MAXCOR)
C
         CALL FSGET(ICORE(I020),1,NUMSYT,1,IRREP,LISTT,'AAII')
         CALL FSGET(ICORE(I010),1,NUMSYL,1,IRREP,LISTL,'IIAI')
C
         CALL XGEMM('N','N',1,NUMSYL,DISSYL,ONEM,ICORE(I020),1,
     $        ICORE(I010),DISSYL,ONE,ICORE(I000),1)
C
C     SEVENTH TERM
C
C     AA AND BB SPINCASE
C
         LISTT=36-ISPIN
         LISTL=138-ISPIN
C
         IRREP=DIRPRD(IRREPO(3-ISPIN),IRREPV(3-ISPIN))
C
         DISSYL=FSDPDII(IRREP,ISYTYP(1,LISTL))
         DISSYT=FSDPDII(IRREP,ISYTYP(1,LISTT))
         NUMSYT=FSDPDAA(IRREP,ISYTYP(2,LISTT))
         NUMSYL=FSDPDAI(IRREP,ISYTYP(2,LISTL))
         NSIZTA=NTIA(3-ISPIN)
         NSIZTB=NTAI(3-ISPIN)
C
C     I000: L2(CONT)  I010: L2   I020:  T2  I030: T1A   I040: T1B
C
         I020=I010+IINTFP*NUMSYL*DISSYL
         I030=I020+IINTFP*NUMSYT*DISSYT
         I040=I030+IINTFP*NSIZTA
         I050=I040+IINTFP*NSIZTB
         IF(I050.GT.MAXCOR) CALL INSMEM('OSL1AI',I050,MAXCOR)
C
         CALL FSGET(ICORE(I020),1,NUMSYT,1,IRREP,LISTT,'IIAA')
         CALL FSGETT1(ICORE(I030),3-ISPIN,90,'IA',11-ISPIN)
         CALL FSGETT1(ICORE(I040),3-ISPIN,90,'AI',11-ISPIN)
         CALL OSTAU35(ICORE(I020),DISSYT,IRREP,
     $        ICORE(I030),NSIZTA,IRREPO(3-ISPIN),VRTI(1,3-ISPIN),
     $        ICORE(I040),NSIZTB,IRREPV(3-ISPIN),POPI(1,3-ISPIN),
     $        ONE,ONE)
         CALL FSGET(ICORE(I010),1,NUMSYL,1,IRREP,LISTL,'IIAI')
C
         CALL XGEMM('T','N',1,NUMSYL,DISSYL,ONEM,ICORE(I020),DISSYT,
     $        ICORE(I010),DISSYL,ONE,ICORE(I000),1)
C
C     AB SPINCASE
C
         LISTT=40-ISPIN
         LISTL=137+ISPIN
C
         IRREP=DIRPRD(IRREPO(ISPIN),IRREPV(3-ISPIN))
C
         DISSYL=FSDPDII(IRREP,ISYTYP(1,LISTL))
         DISSYT=FSDPDII(IRREP,ISYTYP(1,LISTT))
         NUMSYT=FSDPDAA(IRREP,ISYTYP(2,LISTT))
         NUMSYL=FSDPDAI(IRREP,ISYTYP(2,LISTL))
         NSIZTA=NTIA(ISPIN)
         NSIZTB=NTAI(3-ISPIN)
C
C     I000: L2(CONT)  I010: L2   I020:  T2  I030: T1A   I040: T1B
C
         I020=I010+IINTFP*NUMSYL*DISSYL
         I030=I020+IINTFP*NUMSYT*DISSYT
         I040=I030+IINTFP*NSIZTA
         I050=I040+IINTFP*NSIZTB
         IF(I050.GT.MAXCOR) CALL INSMEM('OSL1AI',I050,MAXCOR)
C
         CALL FSGET(ICORE(I020),1,NUMSYT,1,IRREP,LISTT,'IIAA')
         CALL FSGETT1(ICORE(I030),ISPIN,90,'IA',8+ISPIN)
         CALL FSGETT1(ICORE(I040),3-ISPIN,90,'AI',11-ISPIN)
         CALL OSTAU35(ICORE(I020),DISSYT,IRREP,
     $        ICORE(I030),NSIZTA,IRREPO(ISPIN),VRTI(1,ISPIN),
     $        ICORE(I040),NSIZTB,IRREPV(3-ISPIN),POPI(1,3-ISPIN),
     $        ONE,ONE)
         CALL FSGET(ICORE(I010),1,NUMSYL,1,IRREP,LISTL,'IIAI')
C
         CALL XGEMM('T','N',1,NUMSYL,DISSYL,ONEM,ICORE(I020),DISSYT,
     $        ICORE(I010),DISSYL,ONE,ICORE(I000),1)
C
C     SCALE BY W
C
C     I000: L1(CONT)    I010: L1(TARGET)
C
         CALL FSGETT1(ICORE(I010),ISPIN+2,90,'AI',8+ISPIN)
         CALL SAXPY(NSIZTAR,W,ICORE(I000),1,ICORE(I010),1)
         CALL FSPUTT1(ICORE(I010),ISPIN+2,90,'AI',8+ISPIN)
C
 1000 CONTINUE
      RETURN
      END
