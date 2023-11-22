      SUBROUTINE L1INL2A(ICORE,MAXCOR,IUHF,MBPT4)
C
C THIS ROUTINE CALCULATES THE TERM
C
C - SUM M ( <MB//IJ> L(M,A) + <AM//IJ> L(M,B)  
C
C SYMMETRY ADAPTED 
C 
C IN RHF :
C
C - SUM m <Ij//Am> L(m,b) + TRANSPOSITION               
C
C IN UHF
C
C - SUM M <IJ//AM> L(M,B) - ANTISYMMETRIZATION IN A,B
C
C - SUM m <ij//am> L(m,a) - ANTISYMMETRIZATION IN a,b
C
C - SUM m <Ij//Am> L(m,b) - SUM M <Ij//Mb> L(M,A)
C
CEND
C
C CODED AUGUST/90   JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL MBPT4
      INTEGER DIRPRD,DISSYWA,DISSYWB,DISSYZ,POP,VRT,DISSYW
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWO
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
C
      NLIST=100
      IF(MBPT4) NLIST=0
      MXCOR=MAXCOR
      I0TA=MXCOR+1-NTAA*IINTFP
      MXCOR=MXCOR-NTAA*IINTFP
      CALL GETLST(ICORE(I0TA),1,1,1,1,90+NLIST)
      IF(IUHF.EQ.0) THEN
       I0TB=I0TA
      ELSE
       I0TB=I0TA-NTBB*IINTFP
       MXCOR=MXCOR-NTBB*IINTFP
       CALL GETLST(ICORE(I0TB),1,1,1,2,90+NLIST)
      ENDIF
C
C     AA AND BB SPIN CASES
C
      IF(IUHF.EQ.1) THEN
C
C
       DO 100 ISPIN=1,2
C
        IF(ISPIN.EQ.1) THEN
         I000=I0TA
        ELSE
         I000=I0TB
        ENDIF
        LISTZ=ISPIN+60
        LISTW=ISPIN+6
C
        DO 50 IRREP=1,NIRREP 
C
C       RETRIEVE T2 AMPLITUDES AND CALCLUATE NEW ONES
C
         NVRTSQ=0
         DO 45 IRREPJ=1,NIRREP
          IRREPI=DIRPRD(IRREP,IRREPJ)
          NVRTSQ=NVRTSQ+VRT(IRREPJ,ISPIN)*VRT(IRREPI,ISPIN)
45       CONTINUE
         DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
         DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
         NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
         NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
         I001=1
         I002=I001+IINTFP*NVRTSQ*NUMSYZ
         I003=I002+IINTFP*MAX(NUMSYW*DISSYW,NUMSYZ*DISSYZ)
         IF(MIN(NUMSYW,NUMSYZ,DISSYW,DISSYZ).NE.0) THEN
          I004=I003+3*IINTFP*MAX(DISSYW,DISSYZ,NUMSYZ,NUMSYW)
          IF(I004.LT.MXCOR) THEN
C
C
C    IN CORE VERSION
C
          CALL T1T2AA1(ICORE(I001),ICORE(I001),ICORE(I002),ICORE(I002),
     &               ICORE(I002),ICORE(I000),POP(1,ISPIN),VRT(1,ISPIN),
     &               DISSYZ,DISSYW,DISSYW,NUMSYW,NVRTSQ,
     &               LISTZ,LISTW,IRREP,ICORE(I003))
          ELSE
           CALL INSMEM('T1T2AA1',I004,MXCOR)
          ENDIF
         ENDIF
50      CONTINUE
100    CONTINUE
      ENDIF
  
C
C      AB SPIN CASE
C
      LISTZ=63
      LISTWA=9
      IF(IUHF.EQ.0) LISTWA=10
      LISTWB=10
C
C   LOOP OVER IRREPS
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
       I002=I001+IINTFP*NUMSYZ*DISSYZ
       I003=I002+IINTFP*MAX(NUMSYWA*DISSYWA,NUMSYWB*DISSYWB,
     &                      NUMSYZ*DISSYZ)
       IF(MIN(NUMSYZ,DISSYZ).NE.0) THEN
        I004=I003+3*IINTFP*MAX(NUMSYZ,NUMSYWA,NUMSYWB,DISSYWA,
     &                       DISSYWB,DISSYZ)
        IF(I004.LT.MXCOR) THEN
C
C       IN CORE ALGORITHM
C
         CALL T1T2AB1(ICORE(I001),ICORE(I001),ICORE(I002),ICORE(I002),
     &               ICORE(I0TA),ICORE(I0TB),POP(1,1),POP(1,2),VRT(1,1),
     &               VRT(1,2),DISSYZ,DISSYWA,DISSYWB,NUMSYZ,NUMSYWA,
     &               NUMSYWB,NTAA,NTBB,LISTZ,LISTWA,LISTWB,IRREP,IUHF,
     &               ICORE(I003))
        ELSE
C
C       OUT CORE ALGORITHM
C
         CALL INSMEM('T1T2Ab1',I004,MXCOR)
        ENDIF
       ENDIF
200   CONTINUE
      RETURN
      END
