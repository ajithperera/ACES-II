      SUBROUTINE T1INW2_R(ICORE,MAXCOR,IUHF,LISTT,LISTTOFF,IOFFLIST)
C
C THIS ROUTINE CALCULATES THE FOLLOWING CONTRIBUTION TO
C THE W(MN,IJ) INTERMEDIATES :
C
C  P(IJ) SUM E T(E,J) <MN//IE>
C
C THIS CONTRIBUTION IS ONLY REQUIRED FOR CCSD METHODS
C
C RHF :
C
C SUM e <Mn//Ie> T(e,j) + TRANSPOSITION
C
C UHF :
C
C SUM E <MN//IE> T(E,J) + ANTISYMMETRIZATION
C
C SUM e <mn//ie> T(e,j) + ANTISYMMETRIZATION
C
C SUM e <Mn//Ie> T(e,j) + SUM E <Mn//Ej> T(E,I)
C
C IN ADDITION FOR THE F(MI) INTERMEDIATES THE TERMS
C
C SUM N E <MN//IE> T(E,N)
C
C RHF : <MN//IE> T(E,N) + <Mn//Ie> T(e,n)
C
C UHF : <MN//IE> T(E,N) + <Mn//Ie> T(e,n)
C
C       <mn//ie> T(e,n) + <mN//iE> T(E,N)
C
C ARE CALCULATED AS CONTRACTION OF THE PREVIOUSlY
C CALCULATED CONTRIBUTIONS TO W(MN,IJ)
C
CEND
C
C  CODED JULY/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYW,DISSYZ,DISSYWA,DISSYWB,POP,VRT
      DIMENSION ICORE(MAXCOR)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,
     &            NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTIN,IINTFP,IALONE,IBITWO
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /CONTROL/ IPRNT,IXXX,IXXX2
C
C   ALLOCATE CORE MEMORY FOR T1 LIST AND FMI
C
      MXCOR=MAXCOR
      NFAA=NF1AA
      NFBB=NF1BB
      I0AA=MXCOR+1-NFAA*IINTFP
      I0TA=I0AA-NTAA*IINTFP
      MXCOR=MXCOR-(NFAA+NTAA)*IINTFP
      CALL GETLST(ICORE(I0AA),1,1,1,1,91+IOFFLIST)
      CALL GETLST(ICORE(I0TA),1,1,1,1+LISTTOFF,LISTT)
      IF(IUHF.EQ.0) THEN
       I0BB=I0AA
       I0TB=I0TA
      ELSE
       I0BB=I0TA-NFBB*IINTFP
       I0TB=I0BB-NTBB*IINTFP
       MXCOR=MXCOR-(NFBB+NTBB)*IINTFP
       CALL GETLST(ICORE(I0BB),1,1,1,2,91+IOFFLIST)
       CALL GETLST(ICORE(I0TB),1,1,1,2+LISTTOFF,LISTT)
      ENDIF
C
C     AAAA AND BBBB SPIN CASES (UHF ONLY)
C
      IF(IUHF.EQ.1) THEN
C
       DO 1000 ISPIN=1,2
C
        IF(ISPIN.EQ.1) THEN
         IT=I0TA
         IFMI=I0AA
         NT=NTAA
         NF=NFAA
        ELSE
         IT=I0TB
         IFMI=I0BB
         NT=NTBB
         NF=NFBB
        ENDIF
        LISTW=ISPIN+6
        LISTZ=50+ISPIN
C
C    LOOP OVER ALL IRREPS
C
        DO 100 IRREP=1,NIRREP
C
         NOCCSQ=0
         DO 90 IRREPJ=1,NIRREP
          NOCCSQ=NOCCSQ+POP(IRREPJ,ISPIN)*
     &                  POP(DIRPRD(IRREP,IRREPJ),ISPIN)
90       CONTINUE
         DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
         NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
         DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
         NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
         I001=1
         I002=I001+IINTFP*MAX(DISSYW*NUMSYW,NUMSYZ*DISSYZ)
         I003=I002+IINTFP*DISSYW*NOCCSQ
         I004=I003+MAX(NUMSYW,DISSYW,NUMSYZ,DISSYZ)
         IF(MIN(NUMSYW,DISSYW).NE.0) THEN
          IF(I004.LT.MXCOR) THEN
C
C   IN CORE VERSION
C
          CALL T1W2AA_R(ICORE(I001),ICORE(I002),ICORE(IT),ICORE(IFMI),
     &                POP(1,ISPIN),VRT(1,ISPIN),DISSYZ,DISSYW,NUMSYZ,
     &                NUMSYW,NOCCSQ,NT,NF,LISTW,LISTZ,IRREP,ICORE(I003))

          ELSE
C
C    OUT CORE ALGORITHM
C
           STOP 'T1W2AA'
C
          ENDIF
         ENDIF
100     CONTINUE
1000   CONTINUE
      ENDIF
C
C   AB SPIN CASE
C
      LISTWA=10
      LISTWB=9
      LISTZ=53
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
       I002=I001+IINTFP*NUMSYZ*DISSYZ
       I003=I002+IINTFP*MAX(NUMSYWA*DISSYWA,NUMSYWB*DISSYWB,NUMSYZ*
     &                      DISSYZ)
       I004=I003+3*IINTFP*MAX(NUMSYZ,NUMSYWA,NUMSYWB,DISSYZ,DISSYWA,
     &                      DISSYWB)
       IF(I004.LT.MXCOR) THEN
C
        CALL T1W2AB_R(ICORE(I001),ICORE(I002),ICORE(I0TA),ICORE(I0TB),
     &              ICORE(I0AA),ICORE(I0BB),POP(1,1),POP(1,2),VRT(1,1),
     &              VRT(1,2),DISSYZ,DISSYWA,DISSYWB,NUMSYZ,NUMSYWA,
     &              NUMSYWB,NTAA,NTBB,NFAA,NFBB,LISTWA,LISTWB,
     &              LISTZ,IRREP,IUHF,ICORE(I003))
       ELSE
        STOP 'T1W2AB'
       ENDIF
200   CONTINUE
C
C  SAVE FMI CONTRIBUTIONS ON DISK
C
      CALL PUTLST(ICORE(I0AA),1,1,1,1,91+IOFFLIST)
      IF(IUHF.EQ.1) CALL PUTLST(ICORE(I0BB),1,1,1,2,91+IOFFLIST)
      RETURN
      END
