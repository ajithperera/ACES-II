

      SUBROUTINE D2ABIC(UAIA,UAIB,SIJA,SIJB,UIAA,UIAB,SABA,
     &                  SABB,ICORE,MXCOR,IRREPX,IPERT,IUHF,
     &                  STERM,ANTI)
C
C  THIS ROUTINE CALCULATES THE DERIVATIVE INTEGRALS d <Ab|Ic> / d chi
C
C   THERE ARE EIGHT DIFFERENT TERMS:
C
C    d <Ab|Ic>           x         x                 x     
C    --------- = <Ab|Ic>    + SUM U   <Mb/Ic> + SUM U   <Am|Ic>
C      d x                     M   MA            m   mb 
C
C                       x                 x
C                + SUM U   <Ab/Im> + SUM U   <Ab|Ec>
C                   n   mc            E   EI
C
CEND
C
C CODED MARCH/91  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYW,DISSYD,POP,VRT
      LOGICAL FIELD,GEOM,ROHF,QRHF,SEMI,STERM,ANTI
      DIMENSION ICORE(MXCOR) 
      DIMENSION UAIA(1),UAIB(1),UIAA(1),UIAB(1)
      DIMENSION SIJA(1),SIJB(1),SABA(1),SABB(1)
C
      COMMON/OFFSETS/IOFFU(8,2),IOFFS1(8,2),IOFFS2(8,2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPAA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/DTRAN/FIELD,GEOM,ROHF,QRHF,SEMI
C
      DATA AZERO,ONE /0.D0,1.D0/
C
C  LOOP OVER ALL IRREPS OF THE DERIVATIVE INTEGRALS
C
      DO 1000 IRREPR=1,NIRREP
C
C  THE IRREP ON THE LEFT SIDE IS THEN GIVEN AS THE DIRECT PRODUCT
C  OF IRREPP AND IRREPR
C
       IRREPL=DIRPRD(IRREPR,IRREPX)
C
       NUMSYD=IRPDPD(IRREPR,ISYTYP(2,29))
       DISSYD=IRPDPD(IRREPL,ISYTYP(1,29))
C
C  ALLOCATE CORE MEMORY
C
       ID=1
       IW=ID+IINTFP*NUMSYD*DISSYD
C
C  ZERO TARGET LIST
C
       CALL ZERO(ICORE(ID),NUMSYD*DISSYD)
c       call checksum('D2ABIC',icore(Id),numsyd*dissyd)
C
C CONSIDER FIRST TERMS WHICH INVOLVE THE INTEGRALS
C
C   <Am|Ic> and <Mb|Ic>
C
C   <Am/Ic> U(MA) + <Am/Ic> U(mb)
C
C DO <Ic|Mb> U(M,A)
C
        LISTW=26
        NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
        DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
C  ALLOCATE CORE MEMORY
C
        ITMP=IW+IINTFP*NUMSYW*DISSYW
        IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,DISSYD,NUMSYD)
C
        IF(IEND.GE.MXCOR) CALL INSMEM('D2ABIC1',IEND,MXCOR)
C
C  READ IN THE <Ic|Mb> as c,I; b,M
C
        CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPR,LISTW)
C
C  (cI;bA) (R,L) <----- (cI,bM) (R,R) (MA)) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),UIAA,NUMSYD,DISSYD,
     &              DISSYW,NUMSYW,VRT(1,2),POP(1,1),VRT(1,1),
     &              IRREPR,IRREPL,IRREPX,IOFFU(1,1),2)
C
C  TRANSPOSE TARGET LIST cI;bA TO cI;Ab
C
        CALL SYMTR1(IRREPL,VRT(1,2),VRT(1,1),NUMSYD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &              ICORE(ITMP+2*IINTFP*NUMSYD))
c       call checksum('D2ABIC',icore(Id),numsyd*dissyd)
C
C  THE IRREP ON THE RIGHT SIDE IS UNCHANGED, SO IRREPR DETERMINES
C  THE LISTS OF ORIGINAL ZEROTH ORDER INTEGRALS
C
C  <Am|Ic> IS THE SAME AS <Ac|Im>, THUS LIST 22
C 
       LISTW=22
       NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
C ALLOCATE CORE MEMORY
C
       ITMP=IW+IINTFP*NUMSYW*DISSYW
       IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
C
       IF(IEND.GE.MXCOR) CALL INSMEM('D2ABIC2',IEND,MXCOR)
C
C  READ IN THE <Am|Ic> AS <Ac|Im> (Order cI;Am) 
C
        CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPR,LISTW)
C
C  (c,I;A,b)(R,L) <---- (c,I;A,m)(R,R) (m,b) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),UIAB,NUMSYD,DISSYD,
     &              DISSYW,NUMSYW,VRT(1,1),POP(1,2),VRT(1,2),
     &              IRREPR,IRREPL,IRREPX,IOFFU(1,2),2)
C
        IF(STERM) THEN
C
         LISTW=29
         NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
         IW2=IW+IINTFP*NUMSYW*DISSYW
         IEND=IW2+IINTFP*NUMSYW*DISSYW
C
         CALL GETLST(ICORE(IW2),1,NUMSYW,1,IRREPR,LISTW)
C
         CALL TRANSP(ICORE(IW2),ICORE(IW),NUMSYW,DISSYW)
C
         NUMTMP=NUMSYW
         NUMSYW=DISSYW
         DISSYW=NUMTMP
C
         CALL SYMTR3(IRREPR,POP(1,1),VRT(1,2),NUMSYD,DISSYD,
     &               ICORE(ID),ICORE(IW2),ICORE(IW2+IINTFP*DISSYD),
     &               ICORE(IW2+2*IINTFP*DISSYD))
C
C (I,c;A,b) (R,L) <---- (I,c,A,e) (R,R) (e,b) (X)
C
         CALL DFINDT(ICORE(ID),ICORE(IW),SABB,NUMSYD,DISSYD,
     &               DISSYW,NUMSYW,VRT(1,1),VRT(1,2),
     &               VRT(1,2),IRREPR,IRREPL,IRREPX,
     &               IOFFS2(1,2),1)
C
         CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),NUMSYD,ICORE(ID),
     &               ICORE(IW2),ICORE(IW2+IINTFP*NUMSYD),
     &               ICORE(IW2+2*IINTFP*NUMSYD))
         CALL SYMTR1(IRREPR,VRT(1,1),VRT(1,2),DISSYW,ICORE(IW),
     &               ICORE(IW2),ICORE(IW2+IINTFP*DISSYW),
     &               ICORE(IW2+2*IINTFP*DISSYW))

C (I,c;b,A) (R,L) <---- (I,c;b,E) (R,R) (E,A) (X)
C
         CALL DFINDT(ICORE(ID),ICORE(IW),SABA,NUMSYD,DISSYD,
     &              DISSYW,NUMSYW,VRT(1,2),VRT(1,1),VRT(1,1),
     &              IRREPR,IRREPL,IRREPX,IOFFS2(1,1),1)
C
C TRANSPOSE TARGET ARRAY (I,c;b,A) --> (I,c;A,b)
C
        CALL SYMTR1(IRREPL,VRT(1,2),VRT(1,1),NUMSYD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &              ICORE(ITMP+2*IINTFP*NUMSYD))
        ENDIF
C
C  TRANSPOSE TARGET LIST
C
        ID2=ID+NUMSYD*DISSYD*IINTFP
        IEND=ID2+IINTFP*NUMSYD*DISSYD
        IF(IEND.GE.MXCOR) CALL INSMEM('D2ABIC3',IEND,MXCOR)
C
        CALL TRANSP(ICORE(ID),ICORE(ID2),DISSYD,NUMSYD)
c YAU : old
c       CALL ICOPY(IINTFP*NUMSYD*DISSYD,ICORE(ID2),1,ICORE(ID),1)
c YAU : new
        CALL DCOPY(NUMSYD*DISSYD,ICORE(ID2),1,ICORE(ID),1)
c YAU : end
        IF(ANTI) CALL VMINUS(ICORE(ID),NUMSYD*DISSYD)
C
        IF(STERM) THEN
C
         NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
         ITMP=IW+IINTFP*NUMSYW*DISSYW
         IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
         IF(IEND.GE.MXCOR) CALL INSMEM('D2ABIC',IEND,MXCOR)
C
         CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
C  (A,b;I,c) (L,R) <----  (A,b,I,e) (L,L) (e,c) (X)
C
         CALL DFINDT(ICORE(ID),ICORE(IW),SABB,DISSYD,NUMSYD,
     &               DISSYW,NUMSYW,POP(1,1),VRT(1,2),VRT(1,2),
     &              IRREPL,IRREPR,IRREPX,IOFFS2(1,2),1)
C
C TRANSPOSE TARGET AND INTEGRAL LIST A,b;I,c --> A,b;c,I
C
         CALL SYMTR1(IRREPL,POP(1,1),VRT(1,2),DISSYW,ICORE(IW),
     &               ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &               ICORE(ITMP+2*IINTFP*DISSYW))
         CALL SYMTR1(IRREPR,POP(1,1),VRT(1,2),DISSYD,ICORE(ID),
     &               ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYD),
     &               ICORE(ITMP+2*IINTFP*DISSYD))
C
C  (A,b;c,I) (L,R) <---- (A,b,c,M) (L,L) (M,I) (X)
C
         CALL DFINDT(ICORE(ID),ICORE(IW),SIJA,DISSYD,NUMSYD,
     &               DISSYW,NUMSYW,VRT(1,2),POP(1,1),POP(1,1),
     &               IRREPL,IRREPR,IRREPX,IOFFS1(1,1),1)
C
        ENDIF
C
C NOW DEAL WITH THE RIGHT HAND SIDE
C
C  NOW DO <Ab|Ec> U(E,I)
C
       LISTW=233
       NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
C  ALLOCATE MEMORY
C 
       ITMP=IW+IINTFP*NUMSYW*DISSYW
       IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW)
C
       IF(IEND.GT.MXCOR) CALL INSMEM('D2ABIC4',IEND,MXCOR)
C
C  READ IN THE <Ab|Ec>
C
        CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
C  TRANSPOSE INTEGRALS FROM Ab,Ec TO Ab,cE
C
        CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),DISSYW,ICORE(IW),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &              ICORE(ITMP+2*IINTFP*DISSYW))
C
C  (Ab,cI) (L,R) <----- (Ab,cE) (L,L) (E,I) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),UAIA,DISSYD,NUMSYD,
     &              DISSYW,NUMSYW,VRT(1,2),VRT(1,1),POP(1,1),
     &              IRREPL,IRREPR,IRREPX,IOFFU(1,1),1)
C
C TRANSPOSE TARGET LIST
C
        CALL SYMTR1(IRREPR,VRT(1,2),POP(1,1),DISSYD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYD),
     &              ICORE(ITMP+2*IINTFP*DISSYD))
c       call checksum('D2ABIC',icore(Id),numsyd*dissyd)
C
C  NEXT TERM    <Ab,Im> U(mc)
C
       LISTW=16
       DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
       NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
C
C ALLOCATE MEMORY
C
       ITMP=IW+IINTFP*NUMSYW*DISSYW
       IEND=ITMP+3*IINTFP*MAX(NUMSYD,DISSYD,NUMSYW,DISSYW)
       IF(IEND.GT.MXCOR) CALL INSMEM('D2ABIC5',IEND,MXCOR)
C
C READ IN <Ab|Im> INTEGRALS
C
       CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
C  (Ab,Ic) (L,R) <----- (Ab,Im) (L,L) (m,a) (X)
C
       CALL DFINDT(ICORE(ID),ICORE(IW),UIAB,DISSYD,NUMSYD,
     &             DISSYW,NUMSYW,POP(1,1),POP(1,2),VRT(1,2),
     &             IRREPL,IRREPR,IRREPX,IOFFU(1,2),2)
C
       call checksum('D2ABIC',icore(Id),numsyd*dissyd)
       CALL PUTLST(ICORE(ID),1,NUMSYD,1,IRREPR,329)
1000  CONTINUE
C
C ALL DONE SO FAR
C
      RETURN
      END