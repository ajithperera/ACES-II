      SUBROUTINE D2ABCI(UAIA,UAIB,SIJA,SIJB,UIAA,UIAB,SABA,
     &                  SABB,ICORE,MXCOR,IRREPX,IPERT,IUHF,
     &                  STERM,ANTI)
C
C  THIS ROUTINE CALCULATES THE DERIVATIVE INTEGRALS d <Ab|Ci> / d chi
C
C   THERE ARE EIGHT DIFFERENT TERMS:
C
C    d <Ab|Ci>          x         x                 x     
C    --------- = <Ab|Ci>    + SUM U   <Mb|CI> + SUM U   <Am|Ci>
C      d x                     M   MA            m   mb 
C
C                       x                 x
C                + SUM U   <Ab|Mi> + SUM U   <Ab|Ce>
C                   M   MC            e   ei 
C
CEND
C
C CODED MARCH/91  JG
C
C OUT-OF-CORE VERSION OCTOBER 92
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYW,DISSYD,POP,VRT,DISREAD,DISMAX,DISLEFT
      LOGICAL FIELD,GEOM,ROHF,QRHF,SEMI,STERM,ANTI
      LOGICAL NOABCD
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
      COMMON/ABCD/NOABCD
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW 
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      CALL TIMER(1)
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
       NUMSYD=IRPDPD(IRREPR,ISYTYP(2,30))
       DISSYD=IRPDPD(IRREPL,ISYTYP(1,30))
C
C DETERMINE FIRST LENGTH OF OTHER INTEGAL AND SCRATCH ARRAYS
C BEFORE DECIDING ABOUT IN-CORE AND OUT-OF-CORE ALGORITHM.
C
       LISTW=21
       NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
       MEM1=IINTFP*NUMSYW*DISSYW+
     &      3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
       LISTW=25
       NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
       MEM2=IINTFP*NUMSYW*DISSYW+
     &      3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
       LISTW=16
       NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
       MEM3=IINTFP*NUMSYW*DISSYW+
     &      3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
       MEMMAX=MAX(MEM1,MEM2,MEM3)
       MEMTAR=MXCOR-MEMMAX
       IF(MEMTAR.GE.IINTFP*NUMSYD*DISSYD) THEN
C
C  IN-CORE-ALGORITHM
C
C  ALLOCATE CORE MEMORY
C
        ID=1
        IW=ID+IINTFP*NUMSYD*DISSYD
        IF(IW.GE.MXCOR) CALL INSMEM('D2ABCI1',IW,MXCOR)
C
C  ZERO TARGET LIST
C
        CALL ZERO(ICORE(ID),NUMSYD*DISSYD)
C
C CONSIDER FIRST TERMS WHICH INVOLVE THE INTEGRALS
C
C   <Mb|Ci> and <Am|Ci>
C
C   <Mb|Ci> U(MATAR) + <Am|Ci> U(mb)
C
C  THE IRREP ON THE RIGHT SIDE IS UNCHANGED, SO IRREPR DETERMINES
C  THE LISTS OF ORIGINAL ZEROTH ORDER INTEGRALS
C
C  <Mb|Ci> IS THE SAME AS <Cb|Mi>, THUS LIST 21 
C 
        LISTW=21
        NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
        DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
C ALLOCATE CORE MEMORY
C
        ITMP=IW+IINTFP*NUMSYW*DISSYW
        IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
C 
        IF(IEND.GE.MXCOR) CALL INSMEM('D2ABCI2',IEND,MXCOR)
C
C  READ IN THE <Mb|Ci> AS <Cb|Mi> (Order C,i; b,M)
C
        CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPR,LISTW)
C
C  (C,i;b,M)(R,L) <---- (C,i;b,M)(R,R) (M,A) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),UIAA,NUMSYD,DISSYD,
     &              DISSYW,NUMSYW,VRT(1,2),POP(1,1),VRT(1,1),
     &              IRREPR,IRREPL,IRREPX,IOFFU(1,1),2)
C
C  TRANSPOSE THE TARGET LIST  C,i;b,A TO C,i;A,b
C
        CALL SYMTR1(IRREPL,VRT(1,2),VRT(1,1),NUMSYD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &              ICORE(ITMP+2*IINTFP*NUMSYD))
C
C  NOW DO <Ci|Am> U(m,b)
C
        LISTW=25
        NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
        DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
C  ALLOCATE CORE MEMORY
C
        ITMP=IW+IINTFP*NUMSYW*DISSYW
        IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,DISSYD,NUMSYD)
C
        IF(IEND.GE.MXCOR) CALL INSMEM('D2ABCI3',IEND,MXCOR)
C
C  READ IN THE <Ci|Am> as Ci;Am
C
        CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPR,LISTW)
C
C  (Ci;Ab) (R,L) <----- (Ci;Am) (R,R) (m,b) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),UIAB,NUMSYD,DISSYD,
     &              DISSYW,NUMSYW,VRT(1,1),POP(1,2),VRT(1,2),
     &              IRREPR,IRREPL,IRREPX,IOFFU(1,2),2)
C
        IF(STERM) THEN
C
C THE ORIGINAL ALGORITHM REQUIRES THREE TIMES THE LENGTH OF
C THE ABCI LIST, CHANGE TO GETTRN
C
         LISTW=30
         NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
         IW=ID+IINTFP*NUMSYD*DISSYD
         ITMP=IW+IINTFP*NUMSYW*DISSYW
         IEND=ITMP+IINTFP*3*MAX(NUMSYW,NUMSYD,DISSYD,DISSYW)
         IF(IEND.LT.MXCOR) THEN
C
C IN-CORE-ALGORITHM
C
          CALL GETTRN(ICORE(IW),ICORE(ITMP),DISSYW,NUMSYW,
     &                1,IRREPR,LISTW)
          NUMTMP=NUMSYW
          NUMSYW=DISSYW
          DISSYW=NUMTMP
C
C  (Ci;Ab) (R,L) <----- (Ci;Ae) (R,R) (eb) (X)
C
          CALL DFINDT(ICORE(ID),ICORE(IW),SABB,NUMSYD,DISSYD,
     &                DISSYW,NUMSYW,VRT(1,1),VRT(1,2),VRT(1,2),
     &                IRREPR,IRREPL,IRREPX,IOFFS2(1,2),1)
C TRANSPOSE TARGET AND INTEGRAL LISTS
C
          CALL SYMTR1(IRREPR,VRT(1,1),VRT(1,2),DISSYW,ICORE(IW),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &                ICORE(ITMP+2*IINTFP*DISSYW))
C
          CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),NUMSYD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &                ICORE(ITMP+2*IINTFP*NUMSYD))
C
C  (Ci;bA) (R,L) <----- (Ci;bE) (R,R) (EA) (X)
C
          CALL DFINDT(ICORE(ID),ICORE(IW),SABA,NUMSYD,DISSYD,
     &                DISSYW,NUMSYW,VRT(1,2),VRT(1,1),VRT(1,1),
     &                IRREPR,IRREPL,IRREPX,IOFFS2(1,1),1)
C
          CALL SYMTR1(IRREPL,VRT(1,2),VRT(1,1),NUMSYD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &                ICORE(ITMP+2*IINTFP*NUMSYD))
C
         ELSE
C
C OUT-OF-CORE ALGORITHM
C
          ITMP=ID+IINTFP*NUMSYD*DISSYD
          IW=ITMP+IINTFP*3*MAX(NUMSYD,NUMSYW,DISSYD,DISSYW)
          ISTART=1
          NLEFT=NUMSYW
          MEM=MXCOR-IW
          NDIS=MEM/(IINTFP*DISSYW)
          IF(NDIS.EQ.0) THEN
           write(*,*) ' @-D2ABCI-F, Out-of-core-algorithm not possible.'
           CALL ERREX
          ENDIF
900       CONTINUE
C
          NREAD=MIN(NLEFT,NDIS)
          NLEFT=NLEFT-NREAD
          IEND=ISTART+NREAD-1
C
          CALL GETTRN3(ICORE(IW),ICORE(ITMP),DISSYW,
     &                 ISTART,IEND,NREAD,1,IRREPR,LISTW)
C
C  (Ci;Ab) (R,L) <----- (Ci;Eb) (R,R) (EA) (X)
C
          CALL DFINDT2(ICORE(ID),ICORE(IW),SABA,ISTART,IEND,
     &                 NUMSYD,DISSYD,NREAD,DISSYW,VRT(1,2),
     &                 VRT(1,1),VRT(1,1),IRREPR,IRREPL,
     &                 IRREPX,IOFFS2(1,1),1)
C
C TRANSPOSE TARGET AND INTEGRAL LISTS
C   
          CALL SYMTR1(IRREPR,VRT(1,1),VRT(1,2),NREAD,ICORE(IW),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NREAD),
     &                ICORE(ITMP+2*IINTFP*NREAD))
C
          CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),NUMSYD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &                ICORE(ITMP+2*IINTFP*NUMSYD))
C
C  (Ci;bA) (R,L) <----- (Ci;eA) (R,R) (eb) (X)
C
          CALL DFINDT2(ICORE(ID),ICORE(IW),SABB,ISTART,IEND,
     &                 NUMSYD,DISSYD,NREAD,DISSYW,VRT(1,1),
     &                 VRT(1,2),VRT(1,2),IRREPR,IRREPL,
     &                 IRREPX,IOFFS2(1,2),1)
C
          CALL SYMTR1(IRREPL,VRT(1,2),VRT(1,1),NUMSYD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &                ICORE(ITMP+2*IINTFP*NUMSYD))
C
          ISTART=IEND+1
          IF(NLEFT.NE.0) GO TO 900
C
         ENDIF
        ENDIF
C
C TRANSPOSE TARGET ARRAY
C
        ID2=ID+NUMSYD*DISSYD*IINTFP
        IEND=ID2+IINTFP*NUMSYD*DISSYD
        IF(IEND.LT.MXCOR) THEN
C
C IN-CORE-TRANSPOSITION
C 
         CALL TRANSP(ICORE(ID),ICORE(ID2),DISSYD,NUMSYD)
c YAU : old
c        CALL ICOPY(IINTFP*NUMSYD*DISSYD,ICORE(ID2),1,ICORE(ID),1)
c YAU : new
         CALL DCOPY(NUMSYD*DISSYD,ICORE(ID2),1,ICORE(ID),1)
c YAU : end
C
        ELSE
C
C OUT-OF-CORE-TRANSPOSITION
C
         ID2=ID+NUMSYD*DISSYD*IINTFP
         IEND=ID2+IINTFP*NUMSYD
         IF(IEND.GT.MXCOR) CALL INSMEM('D2ABCI4',IEND,MXCOR)
C 
         CALL OOCTRN(ICORE(ID),ICORE(ID),DISSYD,NUMSYD,ICORE(ID2))
C
        ENDIF
C
        IF(ANTI) CALL VMINUS(ICORE(ID),NUMSYD*DISSYD)
C
        IF(STERM) THEN
C 
         LISTW=30
         NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
         IW=ID+IINTFP*NUMSYD*DISSYD
         ITMP=IW+IINTFP*NUMSYW*DISSYW
         IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,DISSYD,NUMSYD)
         IF(IEND.LT.MXCOR) THEN
C
C IN-CORE-VERSION (FOR NOABCD=.FALSE., WE DO NOT SUPPORT 
C OUT-OF-CORE)
C
          CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
          IF(.NOT.NOABCD) THEN
C
C ADD THIS PIECE ONLY IN CASES WHERE IT HAS NOT BEEN CALCULATED
C FROM AO INTEGRALS DIRECTLY
C
C  (Ab,Ci) (L,R) <----- (Ab,Cm) (L,L) (m,i) (X)
C
           CALL DFINDT(ICORE(ID),ICORE(IW),SIJB,DISSYD,NUMSYD,
     &                 DISSYW,NUMSYW,VRT(1,1),POP(1,2),POP(1,2),
     &                 IRREPL,IRREPR,IRREPX,IOFFS1(1,2),1)
C
          ENDIF
C
          CALL SYMTR1(IRREPL,VRT(1,1),POP(1,2),DISSYW,ICORE(IW),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW), 
     &                ICORE(ITMP+2*IINTFP*DISSYW))
C
C TRANSPOSE LAST TWO INDICES OF TARGET ARRAY
C
          CALL SYMTR1(IRREPR,VRT(1,1),POP(1,2),DISSYD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYD), 
     &                ICORE(ITMP+2*IINTFP*DISSYD))
C
C  (Ab,iC) (L,R) <----- (Ab,iE) (L,L) (E,A) (X)
C
          CALL DFINDT(ICORE(ID),ICORE(IW),SABA,DISSYD,NUMSYD,
     &                DISSYW,NUMSYW,POP(1,2),VRT(1,1),VRT(1,1),
     &                IRREPL,IRREPR,IRREPX,IOFFS2(1,1),1)
C
         ELSE
C
C OUT-OF-CORE VERSION
C
          IF(.NOT.NOABCD) THEN
           WRITE(*,*) '@-D2ABCI-F, No out-of-core version available.'
           CALL ERREX
          ENDIF
C
          ITMP=ID+IINTFP*NUMSYD*DISSYD
          IW=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,DISSYD,NUMSYD)
          MEM=MXCOR-IW
          NDIS=MEM/(IINTFP*DISSYW)
          IF(NDIS.LE.0) THEN
           write(*,*) ' @-D2ABCI-F, Out-of-core algorithm not possible.'
           CALL ERREX
          ENDIF
          ISTART=1
          NLEFT=NUMSYW
910       NREAD=MIN(NDIS,NLEFT)
          NLEFT=NLEFT-NREAD
          IEND=ISTART+NREAD-1
C
          CALL GETLST(ICORE(IW),ISTART,NREAD,1,IRREPL,LISTW)
C
C  (Ab,Ci) (L,R) <----- (Ab,Ei) (L,L) (E,A) (X)
C
          CALL DFINDT3(ICORE(ID),ICORE(IW),SABA,ISTART,IEND,
     &                 DISSYD,NUMSYD,DISSYW,NREAD,POP(1,2),
     &                 VRT(1,1),VRT(1,1),IRREPL,IRREPR,
     &                 IRREPX,IOFFS2(1,1),1)
C
          ISTART=IEND+1
          IF(NLEFT.NE.0) GO TO 910
C
C TRANSPOSE LAST TWO INDICES OF TARGET ARRAY
C
          CALL SYMTR1(IRREPR,VRT(1,1),POP(1,2),DISSYD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYD), 
     &                ICORE(ITMP+2*IINTFP*DISSYD))
C
         ENDIF
        ELSE
C
C TRANSPOSE LAST TWO INDICES OF TARGET ARRAY
C
         CALL SYMTR1(IRREPR,VRT(1,1),POP(1,2),DISSYD,ICORE(ID),
     &               ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYD), 
     &               ICORE(ITMP+2*IINTFP*DISSYD))
C
        ENDIF
C
C NOW DEAL WITH THE RIGHT HAND SIDE
C
C  FIRST TERM    <Ab|Mi> U(M,C)
C
        LISTW=16
        DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
        NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
C
C ALLOCATE MEMORY
C
        IW=ID+IINTFP*NUMSYD*DISSYD
        ITMP=IW+IINTFP*NUMSYW*DISSYW
        IEND=ITMP+3*IINTFP*MAX(NUMSYD,DISSYD,NUMSYW,DISSYW)
        IF(IEND.GT.MXCOR) CALL INSMEM('D2ABCI5',IEND,MXCOR)
C
C READ IN <Ab|Mi> INTEGRALS
C
        CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
C TRANSPOSE LAST TWO INDICES: Ab;Mi TO Ab,iM
C
        CALL SYMTR1(IRREPL,POP(1,1),POP(1,2),DISSYW,ICORE(IW),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &              ICORE(ITMP+2*IINTFP*DISSYW))
C
C  (Ab,iC) (L,R) <----- (Ab,iM) (L,L) (M,A) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),UIAA,DISSYD,NUMSYD,
     &              DISSYW,NUMSYW,POP(1,2),POP(1,1),VRT(1,1),
     &              IRREPL,IRREPR,IRREPX,IOFFU(1,1),2)
C
C  TRANSPOSE THE TARGET LIST : A,b;i,C TO A,b;C,i
C
        CALL SYMTR1(IRREPR,POP(1,2),VRT(1,1),DISSYD,
     &              ICORE(ID),ICORE(ITMP),ICORE(ITMP+IINTFP
     &              *DISSYD),ICORE(ITMP+2*IINTFP*DISSYD))
C
C  NOW DO <Ab|Ce> U(ei)
C
        IF(.NOT.NOABCD) THEN
         LISTW=233
         NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
C  ALLOCATE MEMORY
C
         IEND=IW+IINTFP*NUMSYW*DISSYW
C
         IF(IEND.GT.MXCOR) THEN
C
C  OUT OF CORE VERSION
C
          DISMAX=(MXCOR-IW-1)/(2*IINTFP*DISSYW)
          IW1=IW
          IW2=IW1+IINTFP*DISSYW*DISMAX
          IF(DISMAX.EQ.0) THEN
           write(*,*)
     &       ' @-D2ABCI-F, no out-of core algorithm possible'
           CALL ERREX
          ENDIF
          DISLEFT=NUMSYW
          IOFFSET=1
10        CONTINUE
          DISREAD=MIN(DISMAX,DISLEFT)
          DISLEFT=DISLEFT-DISREAD
          CALL GETLST(ICORE(IW2),IOFFSET,DISREAD,1,IRREPL,LISTW) 
          CALL TRANSP(ICORE(IW2),ICORE(IW1),DISREAD,DISSYW)
C
          CALL DFINDT4(ICORE(ID),ICORE(IW1),ICORE(IW2),
     &                 UAIB,DISSYD,NUMSYD,
     &                 DISREAD,DISSYW,IOFFSET,VRT(1,1),VRT(1,2),
     &                 POP(1,2),IRREPL,IRREPR,IRREPX,IOFFU(1,2),1)
C
          IOFFSET=IOFFSET+DISREAD
C
          IF(DISLEFT.NE.0) GO TO 10
C         
         ELSE
C
C IN-CORE VERSION
C
C  READ IN THE <Ab|Ce>
C
          CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
C  (Ab|Ci) (L,R) <----- (Ab,Ce) (L,L) (e,i) (X)
C
          CALL DFINDT(ICORE(ID),ICORE(IW),UAIB,DISSYD,NUMSYD,
     &                DISSYW,NUMSYW,VRT(1,1),VRT(1,2),POP(1,2),
     &                IRREPL,IRREPR,IRREPX,IOFFU(1,2),1)
C
         ENDIF
C
        ENDIF
C
C FOR GEOMETRICAL PERTURBATIONS, INCREMENT INTEGRAL DERIVATIVE LIST
C
        IF(STERM) THEN
C
C GENERAL ALGORITHM (IN-CORE AND OUT_OF-CORE)
C
         IW=ID+IINTFP*NUMSYD*DISSYD
C
         IF(DISSYD.NE.0) THEN
          MEM=MXCOR-IW
          NDIS=MEM/(IINTFP*DISSYD)
          IF(NDIS.LE.0) THEN
           write(*,*) ' @-D2ABCI-F, out-of-core algorithm not possible.'
           CALL ERREX
          ENDIF
          NLEFT=NUMSYD
          ISTART=1
920       NREAD=MIN(NDIS,NLEFT)
          NLEFT=NLEFT-NREAD
          CALL GETLST(ICORE(IW),ISTART,NREAD,1,IRREPR,330)
          CALL SAXPY(NREAD*DISSYD,ONE,ICORE(IW),1,
     &               ICORE(ID+(ISTART-1)*DISSYD*IINTFP),1)
          ISTART=ISTART+NREAD
C 
          IF(NLEFT.NE.0) GO TO 920
C
         ENDIF
        ENDIF
C
        CALL PUTLST(ICORE(ID),1,NUMSYD,1,IRREPR,330)
        call checksum('D2ABCI t',icore(Id),numsyd*dissyd)
C
       ELSE
C
C  OUT-OF-CORE-ALGORITHM
C
C  ALLOCATE CORE MEMORY
C
        NDISD=MEMTAR/(IINTFP*DISSYD)
        ISTARTD=1
        NLEFTD=NUMSYD
9000    CONTINUE
        NTREATD=MIN(NDISD,NLEFTD)
        NLEFTD=NLEFTD-NTREATD
        IENDD=ISTARTD+NTREATD-1
C
        ID=1
        IW=ID+IINTFP*NTREATD*DISSYD
        IF(IW.GE.MXCOR) CALL INSMEM('D2ABCI1',IW,MXCOR)
C
C  ZERO TARGET LIST
C
        CALL ZERO(ICORE(ID),NTREATD*DISSYD)
C
C CONSIDER FIRST TERMS WHICH INVOLVE THE INTEGRALS
C
C   <Mb|Ci> and <Am|Ci>
C
C   <Mb|Ci> U(MA) + <Am|Ci> U(mb)
C
C  THE IRREP ON THE RIGHT SIDE IS UNCHANGED, SO IRREPR DETERMINES
C  THE LISTS OF ORIGINAL ZEROTH ORDER INTEGRALS
C
C  <Mb|Ci> IS THE SAME AS <Cb|Mi>, THUS LIST 21 
C 
        LISTW=21
        NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
        DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
C ALLOCATE CORE MEMORY
C
        ITMP=IW+IINTFP*NUMSYW*NTREATD
        IEND=ITMP+3*IINTFP*MAX(NUMSYW,NTREATD,DISSYD,DISSYW)
C 
        IF(IEND.GE.MXCOR) CALL INSMEM('D2ABCI2',IEND,MXCOR)
C
C  READ IN THE <Mb|Ci> AS <Cb|Mi> (Order C,i; b,M)
C
        CALL GETLST3(ICORE(IW),ICORE(ITMP),ISTARTD,IENDD,
     &               NTREATD,1,NUMSYW,1,IRREPR,LISTW)
C
C  (C,i;b,M)(R,L) <---- (C,i;b,M)(R,R) (M,A) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),UIAA,NTREATD,DISSYD,
     &              NTREATD,NUMSYW,VRT(1,2),POP(1,1),VRT(1,1),
     &              IRREPR,IRREPL,IRREPX,IOFFU(1,1),2)
C
C  TRANSPOSE THE TARGET LIST  C,i;b,A TO C,i;A,b
C
        CALL SYMTR1(IRREPL,VRT(1,2),VRT(1,1),NTREATD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NTREATD),
     &              ICORE(ITMP+2*IINTFP*NTREATD))
C
C  NOW DO <Ci|Am> U(m,b)
C
        LISTW=25
        NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
        DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
C  ALLOCATE CORE MEMORY
C
        ITMP=IW+IINTFP*NUMSYW*NTREATD
        IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,NTREATD,DISSYD)
C
        IF(IEND.GE.MXCOR) CALL INSMEM('D2ABCI3',IEND,MXCOR)
C
C  READ IN THE <Ci|Am> as Ci;Am
C
        CALL GETLST3(ICORE(IW),ICORE(ITMP),ISTARTD,IENDD,NTREATD,
     &               1,NUMSYW,1,IRREPR,LISTW)
C
C  (Ci;Ab) (R,L) <----- (Ci;Am) (R,R) (m,b) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),UIAB,NTREATD,DISSYD,
     &              NTREATD,NUMSYW,VRT(1,1),POP(1,2),VRT(1,2),
     &              IRREPR,IRREPL,IRREPX,IOFFU(1,2),2)
C
        IF(STERM) THEN
C
C THE ORIGINAL ALGORITHM REQUIRES THREE TIMES THE LENGTH OF
C THE ABCI LIST, CHANGE TO GETTRN
C
         LISTW=30
         NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
         IW=ID+IINTFP*NTREATD*DISSYD
         ITMP=IW+IINTFP*NTREATD*DISSYW
         IEND=ITMP+IINTFP*3*MAX(NUMSYW,NTREATD,DISSYD,DISSYW)
         IF(IEND.LT.MXCOR) THEN
C
C IN-CORE-ALGORITHM
C
          CALL GETTRN3(ICORE(IW),ICORE(ITMP),DISSYW,ISTARTD,IENDD,
     &                 NTREATD,1,IRREPR,LISTW)
C
C  (Ci;Ab) (R,L) <----- (Ci;Ae) (R,R) (eb) (X)
C
          CALL DFINDT(ICORE(ID),ICORE(IW),SABB,NTREATD,DISSYD,
     &                NTREATD,DISSYW,VRT(1,1),VRT(1,2),VRT(1,2),
     &                IRREPR,IRREPL,IRREPX,IOFFS2(1,2),1)
C 
C TRANSPOSE TARGET AND INTEGRAL LISTS
C
          CALL SYMTR1(IRREPR,VRT(1,1),VRT(1,2),NTREATD,ICORE(IW),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NTREATD),
     &                ICORE(ITMP+2*IINTFP*NTREATD))
C
          CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),NTREATD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NTREATD),
     &                ICORE(ITMP+2*IINTFP*NTREATD))
C
C  (Ci;bA) (R,L) <----- (Ci;bE) (R,R) (EA) (X)
C
          CALL DFINDT(ICORE(ID),ICORE(IW),SABA,NTREATD,DISSYD,
     &                NTREATD,DISSYW,VRT(1,2),VRT(1,1),VRT(1,1),
     &                IRREPR,IRREPL,IRREPX,IOFFS2(1,1),1)
C
          CALL SYMTR1(IRREPL,VRT(1,2),VRT(1,1),NTREATD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NTREATD),
     &                ICORE(ITMP+2*IINTFP*NTREATD))
C
         ELSE
C
C OUT-OF-CORE ALGORITHM
C
          ITMP=ID+IINTFP*NTREATD*DISSYD
          IW=ITMP+IINTFP*3*MAX(NTREATD,NUMSYW,DISSYD,DISSYW)
          ISTART=ISTARTD
          NLEFT=NTREATD
          MEM=MXCOR-IW
          NDIS=MEM/(IINTFP*DISSYW)
          IF(NDIS.EQ.0) THEN
           write(*,*) ' @-D2ABCI-F, Out-of-core-algorithm not possible.'
           CALL ERREX
          ENDIF
1900      CONTINUE
C
          NREAD=MIN(NLEFT,NDIS)
          NLEFT=NLEFT-NREAD
          IEND=ISTART+NREAD-1
          ISTART1=ISTART-ISTARTD+1
          IEND1=IEND-ISTARTD+1
C
          CALL GETTRN3(ICORE(IW),ICORE(ITMP),DISSYW,
     &                 ISTART,IEND,NREAD,1,IRREPR,LISTW)
C
C  (Ci;Ab) (R,L) <----- (Ci;Eb) (R,R) (EA) (X)
C
          CALL DFINDT2(ICORE(ID),ICORE(IW),SABA,ISTART1,IEND1,
     &                 NTREATD,DISSYD,NREAD,DISSYW,VRT(1,2),
     &                 VRT(1,1),VRT(1,1),IRREPR,IRREPL,
     &                 IRREPX,IOFFS2(1,1),1)
C
C TRANSPOSE TARGET AND INTEGRAL LISTS
C   
          CALL SYMTR1(IRREPR,VRT(1,1),VRT(1,2),NREAD,ICORE(IW),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NREAD),
     &                ICORE(ITMP+2*IINTFP*NREAD))
C
          CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),NTREATD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NTREATD),
     &                ICORE(ITMP+2*IINTFP*NTREATD))
C
C  (Ci;bA) (R,L) <----- (Ci;eA) (R,R) (eb) (X)
C
          CALL DFINDT2(ICORE(ID),ICORE(IW),SABB,ISTART1,IEND1,
     &                 NTREATD,DISSYD,NREAD,DISSYW,VRT(1,1),
     &                 VRT(1,2),VRT(1,2),IRREPR,IRREPL,
     &                 IRREPX,IOFFS2(1,2),1)
C
          CALL SYMTR1(IRREPL,VRT(1,2),VRT(1,1),NTREATD,ICORE(ID),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*NTREATD),
     &                ICORE(ITMP+2*IINTFP*NTREATD))
C
          ISTART=IEND+1
          IF(NLEFT.NE.0) GO TO 1900
C
         ENDIF
        ENDIF
C
C TRANSPOSE TARGET ARRAY
C
        ID2=ID+NTREATD*DISSYD*IINTFP
        IEND=ID2+IINTFP*NTREATD*DISSYD
        IF(IEND.LT.MXCOR) THEN
C
C IN-CORE-TRANSPOSITION
C 
         CALL TRANSP(ICORE(ID),ICORE(ID2),DISSYD,NTREATD)
c YAU : old
c        CALL ICOPY(IINTFP*NTREATD*DISSYD,ICORE(ID2),1,ICORE(ID),1)
c YAU : new
         CALL DCOPY(NTREATD*DISSYD,ICORE(ID2),1,ICORE(ID),1)
c YAU : end
C
        ELSE
C
C OUT-OF-CORE-TRANSPOSITION
C
         ID2=ID+NTREATD*DISSYD*IINTFP
         IEND=ID2+IINTFP*NTREATD
         IF(IEND.GT.MXCOR) CALL INSMEM('D2ABCI4',IEND,MXCOR)
C 
         CALL OOCTRN(ICORE(ID),ICORE(ID),DISSYD,NTREATD,ICORE(ID2))
C
        ENDIF
C
        IF(ANTI) CALL VMINUS(ICORE(ID),NTREATD*DISSYD)
C
        IF(STERM) THEN
C 
         LISTW=30
         NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
         IW=ID+IINTFP*NTREATD*DISSYD
         ITMP=IW+IINTFP*NUMSYW*DISSYW
         IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,DISSYD,NUMSYD)
         IF(IEND.LT.MXCOR) THEN
C
C IN-CORE-VERSION (FOR NOABCD=.FALSE., WE DO NOT SUPPORT 
C OUT-OF-CORE)
C
          CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
          IF(.NOT.NOABCD) THEN
C
C ADD THIS PIECE ONLY IN CASES WHERE IT HAS NOT BEEN CALCULATED
C FROM AO INTEGRALS DIRECTLY
C
C  (Ab,Ci) (L,R) <----- (Ab,Cm) (L,L) (m,i) (X)
C
           WRITE(*,*) '@-D2ABCI-F, No out-of-core version available.'
           CALL ERREX
C
          ENDIF
C
C  (Ab,Ci) (L,R) <----- (Ab,Ei) (L,L) (E,A) (X)
C
          CALL DFINDT6(ICORE(ID),ICORE(IW),SABA,ISTARTD,IENDD,
     &                 DISSYD,NTREATD,DISSYW,NUMSYW,POP(1,2),
     &                 VRT(1,1),VRT(1,1),IRREPL,IRREPR,IRREPX,
     &                 IOFFS2(1,1),1)
C
         ELSE
C
C OUT-OF-CORE VERSION
C
          IF(.NOT.NOABCD) THEN
           WRITE(*,*) '@-D2ABCI-F, No out-of-core version available.'
           CALL ERREX
          ENDIF
C
          ITMP=ID+IINTFP*NTREATD*DISSYD
          IW=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,DISSYD,NTREATD)
          MEM=MXCOR-IW
          NDIS=MEM/(IINTFP*DISSYW)
          IF(NDIS.LE.0) THEN
           write(*,*) ' @-D2ABCI-F, Out-of-core algorithm not possible.'
           CALL ERREX
          ENDIF
          ISTART=1
          NLEFT=NUMSYW
1910      NREAD=MIN(NDIS,NLEFT)
          NLEFT=NLEFT-NREAD
          IEND=ISTART+NREAD-1
C
          CALL GETLST(ICORE(IW),ISTART,NREAD,1,IRREPL,LISTW)
C
C  (Ab,Ci) (L,R) <----- (Ab,Ei) (L,L) (E,A) (X)
C
          CALL DFINDT7(ICORE(ID),ICORE(IW),SABA,ISTARTD,IENDD,
     &                 ISTART,IEND,DISSYD,NTREATD,DISSYW,
     &                 NREAD,POP(1,2),VRT(1,1),VRT(1,1),
     &                 IRREPL,IRREPR,IRREPX,IOFFS2(1,1),1)
C
          ISTART=IEND+1
          IF(NLEFT.NE.0) GO TO 1910
C
          ENDIF
         ELSE
C
         ENDIF
C
C NOW DEAL WITH THE RIGHT HAND SIDE
C
C  FIRST TERM    <Ab|Mi> U(M,C)
C
         LISTW=16
         DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
         NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
C
C ALLOCATE MEMORY
C
         IW=ID+IINTFP*NTREATD*DISSYD
         ITMP=IW+IINTFP*NUMSYW*DISSYW
         IEND=ITMP+3*IINTFP*MAX(NTREATD,DISSYD,NUMSYW,DISSYW)
         IF(IEND.GT.MXCOR) CALL INSMEM('D2ABCI5',IEND,MXCOR)
C
C READ IN <Ab|Mi> INTEGRALS
C
         CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
C  (Ab,Ci) (L,R) <----- (Ab,Mi) (L,L) (M,A) (X)
C
         CALL DFINDT6(ICORE(ID),ICORE(IW),UIAA,ISTARTD,IENDD,
     &                DISSYD,NTREATD,DISSYW,NUMSYW,POP(1,2),
     &                POP(1,1),VRT(1,1),
     &                IRREPL,IRREPR,IRREPX,IOFFU(1,1),2)
C
C  NOW DO <Ab|Ce> U(ei)
C
         IF(.NOT.NOABCD) THEN
C
C  OUT OF CORE VERSION
C
          write(*,*)
     &      ' @-D2ABCI-F, no out-of core algorithm possible'
          CALL ERREX
C
         ENDIF
C
C FOR GEOMETRICAL PERTURBATIONS, INCREMENT INTEGRAL DERIVATIVE LIST
C
         IF(STERM) THEN
C
C GENERAL ALGORITHM (IN-CORE AND OUT_OF-CORE)
C
          IW=ID+IINTFP*DISSYD*NTREATD
C
          IF(DISSYD.NE.0) THEN
           MEM=MXCOR-IW
           NDIS=MEM/(IINTFP*DISSYD)
           IF(NDIS.LE.0) THEN
            write(*,*) ' @-D2ABCI-F, out-of-core algorithm',
     &                 ' not possible.'
            CALL ERREX
           ENDIF
           NLEFT=NTREATD
           ISTART=ISTARTD
1920       NREAD=MIN(NDIS,NLEFT)
           NLEFT=NLEFT-NREAD
           CALL GETLST(ICORE(IW),ISTART,NREAD,1,IRREPR,330)
           CALL SAXPY(NREAD*DISSYD,ONE,ICORE(IW),1,
     &                ICORE(ID+(ISTART-ISTARTD)*DISSYD*IINTFP),1)
           ISTART=ISTART+NREAD
C 
           IF(NLEFT.NE.0) GO TO 1920
C
          ENDIF
         ENDIF
C
         CALL PUTLST(ICORE(ID),ISTARTD,NTREATD,1,IRREPR,330)
         call checksum('D2ABCI t',icore(Id),ntreatd*dissyd)
C
         ISTARTD=IENDD+1
         IF(NLEFTD.NE.0) GO TO 9000
C
         ENDIF
      
1000    CONTINUE
C
C ALL DONE SO FAR, WRITE INFO MESSAGE AND GET CPU TIMING
C
       CALL TIMER(1) 
       write(6,6000) TIMENEW
6000   FORMAT(' Integral derivatives d <ab||ci>/d chi',
     &            ' have been formed in ',f5.1,' seconds.')
C
      RETURN
      END
