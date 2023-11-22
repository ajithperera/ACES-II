      SUBROUTINE FORMW1(ICORE,MAXCOR,IUHF,ADD)
C
C THIS ROUTINE CALCULATES THE FOLLOWING CONTRIBUTION TO 
C THE W(AB,EF) INTERMEDIATES :
C
C - P(AB) SUM M T(M,B) <AM//EF>
C
C THIS CONTRIBUTION IS ONLY REQUIRED FOR CCSD METHODS
C
C RHF :
C
C - SUM m <Ef//Am> T(m,b) - TRANSPOSITION
C
C UHF :
C
C - SUM M <EF//AM> T(M,B) - ANTISYMMETRIZATION
C
C - SUM m <ef//am> T(m,b) - ANTISYMMETRIZATION
C
C - SUM m <Ef//Am> T(m,b) - SUM M <Ef//Mb> T(M,A)
C
CEND
C
C  CODED SEPTEMBER/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL ADD
      INTEGER DIRPRD,DISSYW,DISSYZ,DISSYWA,DISSYWB,POP,VRT
      DIMENSION ICORE(MAXCOR),I0T(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /MACHSP/ IINTLN,IFLTIN,IINTFP,IALONE,IBITWO
C
      DATA ONE,ONEM /1.0D0,-1.D0/
C
      FACT=ONE
      IF(ADD) FACT=ONEM
C
C   ALLOCATE CORE MEMORY FOR T1 LIST AND FMI
C
      MXCOR=MAXCOR
      I0T(1)=MXCOR+1-NT(1)*IINTFP
      MXCOR=MXCOR-NT(1)*IINTFP
      CALL GETLST(ICORE(I0T(1)),1,1,1,1,90)
      IF(IUHF.EQ.0) THEN
       I0T(2)=I0T(1)
      ELSE
       I0T(2)=I0T(1)-NT(2)*IINTFP
       MXCOR=MXCOR-NT(2)*IINTFP
       CALL GETLST(ICORE(I0T(2)),1,1,1,2,90)
      ENDIF 
C
C     AAAA AND BBBB SPIN CASES (UHF ONLY)
C
      IF(IUHF.EQ.1) THEN
C
       DO 1000 ISPIN=1,2
C
        LISTW=ISPIN+26
        LISTZ=230+ISPIN
C
C    LOOP OVER ALL IRREPS
C
        DO 100 IRREP=1,NIRREP
C
         NVRTSQ=0
         DO 90 IRREPJ=1,NIRREP
          NVRTSQ=NVRTSQ+VRT(IRREPJ,ISPIN)*
     &                  VRT(DIRPRD(IRREP,IRREPJ),ISPIN)
90       CONTINUE
         DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
         NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
         DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
         NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
         I001=1
         I002=I001+IINTFP*DISSYW*NUMSYW
         I003=I002+MAX(NUMSYW,DISSYW,NUMSYZ,DISSYZ)*IINTFP
         MAXSIZE=(MXCOR-I003)/IINTFP
         IF(MIN(NUMSYW,DISSYW).NE.0) THEN
          IF(MAXSIZE.GT.DISSYZ) THEN
C
C   IN CORE VERSION
C
          CALL W1AA(ICORE(I001),ICORE(I003),MAXSIZE,
     &              ICORE(I0T(ISPIN)),FACT,VRT(1,ISPIN),
     &              POP(1,ISPIN),DISSYZ,DISSYW,NUMSYZ,NUMSYW,
     &              NVRTSQ,LISTW,LISTZ,IRREP,ICORE(I002))
                    
          ELSE
C
C    OUT CORE ALGORITHM
C
           WRITE(*,*)'STOP W1AA'
           CALL ERREX
C
          ENDIF
         ENDIF
100     CONTINUE
1000   CONTINUE
      ENDIF
C
C   AB SPIN CASE
C
      LISTWA=30
      LISTWB=29
      LISTZ=233
C
C     LOOP OVER IRREPS
C
      DO 200 IRREP=1,NIRREP
C
       DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
       DISSYWA=IRPDPD(IRREP,ISYTYP(1,LISTWA))
       NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
       NUMSYWA=IRPDPD(IRREP,ISYTYP(2,LISTWA))
       IF(IUHF.EQ.0) THEN
        DISSYWB=0
        NUMSYWB=0
       ELSE
        DISSYWB=IRPDPD(IRREP,ISYTYP(1,LISTWB))
        NUMSYWB=IRPDPD(IRREP,ISYTYP(2,LISTWB))
       ENDIF
       I001=1
       I002=I001+IINTFP*MAX(NUMSYWA*DISSYWA,NUMSYWB*DISSYWB)
       I003=I002+3*IINTFP*MAX(NUMSYZ,NUMSYWA,NUMSYWB,DISSYZ,DISSYWA,
     &                      DISSYWB)
       MAXSIZE=(MXCOR-I003)/IINTFP
       IF(MAXSIZE.GT.DISSYWA) THEN
C
        CALL W1AB(ICORE(I001),ICORE(I003),MAXSIZE,ICORE(I0T(1)),
     &            ICORE(I0T(2)),FACT,ADD,VRT(1,1),VRT(1,2),POP(1,1),
     &            POP(1,2),DISSYZ,DISSYWA,DISSYWB,NUMSYZ,NUMSYWA,
     &            NUMSYWB,LISTWA,LISTWB,LISTZ,IRREP,IUHF,ICORE(I002))
       ELSE
        WRITE(*,*) 'STOP W1AB'
        CALL ERREX
       ENDIF
200   CONTINUE
      RETURN
      END
