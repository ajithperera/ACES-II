 
      SUBROUTINE H4G5AB(H1,H2,G,G1,T1A,T1B,FACT,DISSYH1,DISSYH2,
     &                  DISSYG,NUMSYH1,NUMSYH2,NUMSYG,POP1,POP2,    
     &                  VRT1,VRT2,LISTH1,LISTH2,LISTG,ISPIN,IRREP,
     &                  IUHF,TMP,TRIP,LISTG1,INCORE)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL TRIP,INCORE,CIS,EOM,LTRP
      INTEGER DISSYH1,DISSYH2,DISSYG,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION H1(DISSYH1,NUMSYH1),H2(DISSYH2,NUMSYH2),
     &          G(NUMSYG,DISSYG),G1(DISSYG,*),T1A(1),T1B(1),
     &          TMP(1),POP1(8),POP2(8),VRT1(8),VRT2(8) 
C
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/EXCITE/CIS,EOM
      COMMON/LTRIP/LTRP
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      CALL ZERO(G,NUMSYG*DISSYG)
C
C  GET RESORTED H4 FROM DISK: FIRST THE H4(aN,iC) INTERMEDIATES
C
      CALL GETLST(H1,1,NUMSYH1,1,IRREP,LISTH1)
C
C  THE ORDER IS ALREADY: Ci, bN

C  LOOP OVER IRREPS OF N AND MULTIPLY WITH T1A
C
      IOFFH=1
      IOFFG=1
      IOFFT=1
      DO 100 IRREPJ=1,NIRREP
C
C  GET POPULATIONS
C
       NVRTJ=VRT1(IRREPJ)
       NOCCJ=POP1(IRREPJ)
       IRREPI=DIRPRD(IRREPJ,IRREP)
       NVRTI=VRT2(IRREPI)
C
C IF ONE OF THE POPULATIONS IS ZERO, SKIP MULTIPLICATION 
C
       IF(MIN(NVRTJ,NVRTI,NOCCJ).NE.0) THEN
C
        CALL XGEMM('N','T',DISSYH1*NVRTI,NVRTJ,NOCCJ,FACT,H1(1,IOFFH),
     &             DISSYH1*NVRTI,T1A(IOFFT),NVRTJ,AZERO,G(1,IOFFG),
     &             NUMSYG*NVRTI)
       ENDIF
C
C UPDATE OFFSETS
C
       IOFFT=IOFFT+NOCCJ*NVRTJ
       IOFFG=IOFFG+NVRTI*NVRTJ
       IOFFH=IOFFH+NVRTI*NOCCJ
C
100    CONTINUE 
C
C  TRANSPOSE THE G LIST FOR SECOND CONTRACTION :  C,i ; b,A --> C,I ; A,b
C
       CALL SYMTR1(IRREP,VRT2,VRT1,NUMSYG,G,TMP,TMP(1+NUMSYG),
     &             TMP(1+2*NUMSYG))
C
C  GET SECOND H4 LIST : H4(An,Ci) INTERMEDIATE WITH ORDERING A,n ; C,i
C
       CALL GETLST(H2,1,NUMSYH2,1,IRREP,LISTH2)
C
C   THE ORDER IS ALREADY Ci,An
C
C      PERFORM MULTIPLICATION WITH T1B
C
      IOFFG=1
      IOFFH=1
      IOFFT=1
      DO 200 IRREPJ=1,NIRREP
C
C  GET POPULATIONS FOR MULTIPLICATION
C
       NOCCJ=POP2(IRREPJ)
       NVRTJ=VRT2(IRREPJ)
       IRREPI=DIRPRD(IRREP,IRREPJ)
       NVRTI=VRT1(IRREPI)
C
C  IF ONE OF THE POPULATIONS IS ZERO, SKIP THE MULTIPLICATION 
C
       IF(MIN(NVRTI,NOCCJ,NVRTJ).NE.0) THEN
C
        CALL XGEMM('N','T',DISSYH2*NVRTI,NVRTJ,NOCCJ,FACT,H2(1,IOFFH),
     &             DISSYH2*NVRTI,T1B(IOFFT),NVRTJ,ONE,G(1,IOFFG),
     &             NUMSYG*NVRTI)
C
       ENDIF
C
C   UPDATE OFFSETS
C
       IOFFT=IOFFT+NOCCJ*NVRTJ
       IOFFG=IOFFG+NVRTI*NVRTJ
       IOFFH=IOFFH+NVRTI*NOCCJ
C
200   CONTINUE
C
       IF(ISPIN.EQ.2) THEN
        CALL SYMTR1(IRREP,VRT1,VRT2,NUMSYG,G,TMP,TMP(1+NUMSYG),
     &              TMP(1+2*NUMSYG))
       ENDIF
C
C  TRANPOSE FIRST TWO INDICES G(cI,Ab)->G(cI,Ab)
C
       IF(ISPIN.EQ.2) THEN
        CALL SYMTR3(IRREP,VRT1,POP2,NUMSYG,DISSYG,G,TMP,TMP(1+DISSYG),
     &              TMP(1+2*DISSYG))
       ENDIF
C
       IF(INCORE)THEN
C
C  TRANSPOSE THE WHOLE G MATRIX : C,i ; A,b
C
        CALL TRANSP(G,G1,DISSYG,NUMSYG)
C
C IN THE CASE OF TRIPLE EXCITATIONS ADD HERE THEIR CONTRIBUTION
C
        IF(TRIP) THEN
         CALL GETLST(G,1,NUMSYG,1,IRREP,LISTG1)
         CALL SAXPY(NUMSYG*DISSYG,ONE,G,1,G1,1)
        ENDIF
        IF(EOM .OR. LTRP) THEN
         CALL GETLST(G,1,NUMSYG,1,IRREP,LISTG)
         CALL SAXPY(NUMSYG*DISSYG,ONE,G,1,G1,1)
        ENDIF
C 
C SAVE GAMMA6 ON DISK
C
        CALL PUTLST(G1,1,NUMSYG,1,IRREP,LISTG)
C
       ELSE
C
        DO 300 IDIS=1,NUMSYG
         CALL SCOPY(DISSYG,G(IDIS,1),NUMSYG,G1,1)
         IF(TRIP)THEN
          CALL GETLST(G1(1,2),IDIS,1,1,IRREP,LISTG1)
          CALL SAXPY (DISSYG,ONE,G1(1,2),1,G1,1)
         ENDIF
         IF(EOM .OR. LTRP)THEN
          CALL GETLST(G1(1,2),IDIS,1,1,IRREP,LISTG)
          CALL SAXPY (DISSYG,ONE,G1(1,2),1,G1,1)
         ENDIF
         CALL PUTLST(G1,IDIS,1,1,IRREP,LISTG)
300     CONTINUE
C
       ENDIF        
C   
C  ALL DONE, RETURN
C
       RETURN
       END
