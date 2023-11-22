
      SUBROUTINE D2IJKA(UAIA,UAIB,SIJA,SIJB,UIAA,UIAB,SABA,
     &                  SABB,ICORE,MXCOR,IRREPX,IPERT,IUHF,
     &                  STERM,ANTI)
C
C  THIS ROUTINE CALCULATES THE DERIVATIVE INTEGRALS d <Ij|Ka> / d chi
C
C   THERE ARE EIGHT DIFFERENT TERMS:
C
C    d <Ij/Ka>           x         x                 x     
C    --------- = <Ij/Ka>    + SUM U   <Ej/Ka> + SUM U   <Ie/Ka>
C      d x                     E   EI            e   ej 
C
C                       x                 x
C                + SUM U   <Ij/Ea> + SUM U   <Ij/Km>
C                   E   EK            m   ma
C
C
C  FOR UHF THERE TWO SPIN CASES, WHILE FOR RHF THERE IS ONLY ONE
C
C
CEND
C
C CODED MARCH/91  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYW,DISSYD,POP,VRT
      LOGICAL STERM,ANTI
      DIMENSION ICORE(MXCOR) 
      DIMENSION UAIA(1),UAIB(1),UIAA(1),UIAB(1)
      DIMENSION SIJA(1),SIJB(1),SABA(1),SABB(1)
C
      COMMON/OFFSETS/IOFFU(8,2),IOFFS1(8,2),IOFFS2(8,2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPAA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW 
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      CALL TIMER(1)
C
C  INCORE ALGORITHM
C
C  DO HERE ALPHA,BETA SPIN CASE
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
       NUMSYD=IRPDPD(IRREPR,ISYTYP(2,10))
       DISSYD=IRPDPD(IRREPL,ISYTYP(1,10))
C
C  ALLOCATE CORE MEMORY
C
       ID=1
       IW=ID+IINTFP*NUMSYD*DISSYD
       IF(IW.GE.MXCOR) CALL INSMEM('D2IJKA1',IW,MXCOR)
C
C  ZERO TARGET LIST
C
       CALL ZERO(ICORE(ID),NUMSYD*DISSYD)
C
C CONSIDER FIRST TERMS WHICH INVOLVE THE INTEGRALS
C
C   <Ej/Ka> and <Ie/Ka>
C
C   <Ej/Ka> U(EI) + <Ie/Ka> U(ej)
C
C  THE IRREP ON THE RIGHT SIDE IS UNCHANGED, SO IRREPR DETERMINES
C  THE LISTS OF ORIGINAL ZEROTH ORDER INTEGRALS
C
C  <Ej/Ka> IS THE SAME AS <Ea/Kj>, THUS LIST 21 
C 
       LISTW=21
       NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
C ALLOCATE CORE MEMORY
C
       ITMP=IW+IINTFP*NUMSYW*DISSYW
       IW2=ITMP
       IEND1=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
       IEND2=IW2+IINTFP*NUMSYW*DISSYW
       IEND=MAX(IEND1,IEND2)
C
       IF(IEND.GE.MXCOR) THEN
C
        IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
        IF(IEND.GE.MXCOR) CALL INSMEM('D2IJKA2',IEND,MXCOR)
        CALL GETTRN(ICORE(IW),ICORE(ITMP),DISSYW,NUMSYW,1,IRREPR,LISTW)
C
       ELSE
C
C  READ IN THE <Ej/Ka> AS <Ea/Kj> (Order E,j; a,K)
C
        CALL GETLST(ICORE(IW2),1,NUMSYW,1,IRREPR,LISTW)
C
C  TRANSPOSE FROM E,j ; a,K to a,K ; E,j
C
        CALL TRANSP(ICORE(IW2),ICORE(IW),NUMSYW,DISSYW)
C
       ENDIF
C
       NUMTMP=NUMSYW
       NUMSYW=DISSYW
       DISSYW=NUMTMP
C
C  TRANSPOSE FROM a,K ; E,j TO a,K ; j,E
C
       CALL SYMTR1(IRREPR,VRT(1,1),POP(1,2),DISSYW,ICORE(IW),
     &             ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &             ICORE(ITMP+2*IINTFP*DISSYW))
C
C  (a,K;j,I)(R,L) <---- (a,K;j,E)(R,R) (E,I) (X)
C
       CALL DFINDT(ICORE(ID),ICORE(IW),UAIA,NUMSYD,DISSYD,
     &             DISSYW,NUMSYW,POP(1,2),VRT(1,1),POP(1,1),
     &             IRREPR,IRREPL,IRREPX,IOFFU(1,1),1)
C
C  TRANSPOSE THE TARGET LIST  a,K;j,I TO a,k;I,j
C
       CALL SYMTR1(IRREPL,POP(1,2),POP(1,1),NUMSYD,ICORE(ID),
     &             ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &             ICORE(ITMP+2*IINTFP*NUMSYD))
C
C  NOW DO <Ka|Ie> U(ej)
C
       LISTW=25+IUHF
       NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
C  ALLOCATE CORE MEMORY
C
       ITMP=IW+IINTFP*NUMSYW*DISSYW
       IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,DISSYD,NUMSYD)
C
       IF(IEND.GE.MXCOR) CALL INSMEM('D2IJKA3',IEND,MXCOR)
C
C  READ IN THE <Ka|Ie> as a,K; e,I
C
       CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPR,LISTW)
C
C  TRANSPOSE a,K; e,I TO a,K; I,e
C
       CALL SYMTR1(IRREPR,VRT(1,2),POP(1,1),DISSYW,ICORE(IW),
     &             ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &             ICORE(ITMP+2*IINTFP*DISSYW))
C
C  (aK;Ij) (R,L) <----- (aK,Ie) (R,R) (e,j) (X)
C
       CALL DFINDT(ICORE(ID),ICORE(IW),UAIB,NUMSYD,DISSYD,
     &             DISSYW,NUMSYW,POP(1,1),VRT(1,2),POP(1,2),
     &             IRREPR,IRREPL,IRREPX,IOFFU(1,2),1)
C
C  TRANSPOSE TARGET LIST
C
C  ADD SOME ADDITIONAL S-TERMS
C
       IF(STERM) THEN
C
C    FIRST    <Ka|Im> S(mj) + <Ka|Mj> s(MI)
C
        LISTW=10
        NUMSYW=IRPDPD(IRREPR,ISYTYP(2,LISTW))
        DISSYW=IRPDPD(IRREPR,ISYTYP(1,LISTW))
C
        IW=ID+NUMSYD*DISSYD*IINTFP
        IW2=IW+NUMSYW*DISSYW*IINTFP
        ITMP=IW2
        IEND1=IW2+IINTFP*NUMSYW*DISSYW
        IEND2=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
        IEND=MAX(IEND1,IEND2)
        IF(IEND.GE.MXCOR) THEN
C
         IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
         IF(IEND.GE.MXCOR) CALL INSMEM('D2IJKA4',IEND,MXCOR)
         CALL GETTRN(ICORE(IW),ICORE(ITMP),DISSYW,NUMSYW,1,IRREPR,LISTW)
C
        ELSE
C
         CALL GETLST(ICORE(IW2),1,NUMSYW,1,IRREPR,LISTW)
         CALL TRANSP(ICORE(IW2),ICORE(IW),NUMSYW,DISSYW)
C
        ENDIF
C
        NUMTMP=NUMSYW
        NUMSYW=DISSYW
        DISSYW=NUMTMP

        CALL SYMTR3(IRREPR,POP(1,1),VRT(1,2),DISSYW,NUMSYW,
     &              ICORE(IW),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYW),
     &              ICORE(ITMP+2*IINTFP*NUMSYW))
C
C   (aK,Ij) (RL) <---- (aK,Im) (R,R) S(mj) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),SIJB,NUMSYD,DISSYD,
     &              DISSYW,NUMSYW,POP(1,1),POP(1,2),POP(1,2),
     &              IRREPR,IRREPL,IRREPX,IOFFS1(1,2),1)
C
C TRANSPOSE BOTH TARGET AND INTEGRAL LIST
C
        CALL SYMTR1(IRREPR,POP(1,1),POP(1,2),DISSYW,ICORE(IW),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &              ICORE(ITMP+2*IINTFP*DISSYW))
C
        CALL SYMTR1(IRREPL,POP(1,1),POP(1,2),NUMSYD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &              ICORE(ITMP+2*IINTFP*NUMSYD))
C
C   (aK,jI) (RL) <---- (aK,jM) (R,R) S(MI) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),SIJA,NUMSYD,DISSYD,
     &              DISSYW,NUMSYW,POP(1,2),POP(1,1),POP(1,1),
     &              IRREPR,IRREPL,IRREPX,IOFFS1(1,1),1)
C
C TRANSPOSE TARGET LIST
C
        CALL SYMTR1(IRREPL,POP(1,2),POP(1,1),NUMSYD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*NUMSYD),
     &              ICORE(ITMP+2*IINTFP*NUMSYD))
C
       ENDIF
C        
       ID2=ID+NUMSYD*DISSYD*IINTFP
       IEND=ID2+IINTFP*NUMSYD*DISSYD
       IF(IEND.GE.MXCOR) CALL INSMEM('D2IJKA5',IEND,MXCOR)
C
       CALL TRANSP(ICORE(ID),ICORE(ID2),DISSYD,NUMSYD)
c YAU : old
c      CALL ICOPY(IINTFP*NUMSYD*DISSYD,ICORE(ID2),1,ICORE(ID),1)
c YAU : new
       CALL DCOPY(NUMSYD*DISSYD,ICORE(ID2),1,ICORE(ID),1)
c YAU : end
       IF(ANTI) CALL VMINUS(ICORE(ID),NUMSYD*DISSYD)
C
       IF(STERM) THEN
C
        LISTW=10
        NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
        DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
        IW=ID+IINTFP*NUMSYD*DISSYD
        ITMP=IW+IINTFP*NUMSYW*DISSYW
        IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,NUMSYD,DISSYD)
        IF(IEND.GE.MXCOR) CALL INSMEM('D2IJKA6',IEND,MXCOR)
C
        CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
        CALL SYMTR1(IRREPL,POP(1,1),VRT(1,2),DISSYW,ICORE(IW),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &              ICORE(ITMP+IINTFP*2*DISSYW))
C
C  (Ij,aK) (L,R) <--- (Ij,aM) S(MK) (L,L) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),SIJA,DISSYD,NUMSYD,
     &              DISSYW,NUMSYW,VRT(1,2),POP(1,1),POP(1,1),
     &              IRREPL,IRREPR,IRREPX,IOFFS1(1,1),1)
C
C  TRANSPOSE BOTH TARGET AND INTEGRAL ARRAY
C
        CALL SYMTR1(IRREPL,VRT(1,2),POP(1,1),DISSYW,ICORE(IW),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &              ICORE(ITMP+IINTFP*2*DISSYW))
C
        CALL SYMTR1(IRREPR,VRT(1,2),POP(1,1),DISSYD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYD),
     &              ICORE(ITMP+IINTFP*2*DISSYD))
C
C  (Ij,Ka) (L,R) <--- (Ij,Ke) S(ea) (L,L) (X)
C
        CALL DFINDT(ICORE(ID),ICORE(IW),SABB,DISSYD,NUMSYD,
     &              DISSYW,NUMSYW,POP(1,1),VRT(1,2),VRT(1,2),
     &              IRREPL,IRREPR,IRREPX,IOFFS2(1,1),1)
C
        CALL SYMTR1(IRREPR,POP(1,1),VRT(1,2),DISSYD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYD),
     &              ICORE(ITMP+IINTFP*2*DISSYD))
C
       ENDIF
C
C NOW DEAL WITH THE RIGHT HAND SIDE
C
C  FIRST TERM    <Ij|Ea> U(E,K)
C
       LISTW=16
       DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
       NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
C
C ALLOCATE MEMORY
C
       ITMP=IW+IINTFP*NUMSYW*DISSYW
       IW2=ITMP
       IEND1=ITMP+3*IINTFP*MAX(NUMSYD,DISSYD,NUMSYW,DISSYW)
       IEND2=IW2+IINTFP*NUMSYW*DISSYW
       IEND=MAX(IEND1,IEND2)
       IF(IEND.GE.MXCOR) THEN
C
        IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,DISSYD,NUMSYD)
        IF(IEND.GE.MXCOR) CALL INSMEM('D2IJKA7',IEND,MXCOR)
        CALL GETTRN(ICORE(IW),ICORE(ITMP),DISSYW,NUMSYW,1,IRREPL,LISTW)
C
       ELSE
C
C READ IN <Ij|Ea> INTEGRALS
C
        CALL GETLST(ICORE(IW2),1,NUMSYW,1,IRREPL,LISTW)
C 
C TRANSPOSE INTEGRALS Ea,Ij TO Ij,Ea
C
        CALL TRANSP(ICORE(IW2),ICORE(IW),NUMSYW,DISSYW)
C
       ENDIF
C
       NUMTMP=NUMSYW
       NUMSYW=DISSYW
       DISSYW=NUMTMP
C
C TRANSPOSE LAST TWO INDICES:  Ij,Ea TO Ij,aE
C
       CALL SYMTR1(IRREPL,VRT(1,1),VRT(1,2),DISSYW,ICORE(IW),
     &             ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &             ICORE(ITMP+2*IINTFP*DISSYW))
C
C  (Ij,aK) (L,R) <----- (Ij,aE) (L,L) (E,I) (X)
C
       CALL DFINDT(ICORE(ID),ICORE(IW),UAIA,DISSYD,NUMSYD,
     &             DISSYW,NUMSYW,VRT(1,2),VRT(1,1),POP(1,1),
     &             IRREPL,IRREPR,IRREPX,IOFFU(1,1),1)
C
C  TRANSPOSE THE TARGET LIST : Ij,aK TO Ij,Ka
C
       CALL SYMTR1(IRREPR,VRT(1,2),POP(1,1),DISSYD,
     &             ICORE(ID),ICORE(ITMP),ICORE(ITMP+IINTFP
     &             *DISSYD),ICORE(ITMP+2*IINTFP*DISSYD))
C
C  NOW DO <Ij|Km> U(ma)
C
       LISTW=13
       NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
C  ALLOCATE MEMORY
C
       IEND=IW+IINTFP*NUMSYW*DISSYW
C
       IF(IEND.GT.MXCOR) CALL INSMEM('D2IJKA8',IEND,MXCOR)
C
C  READ IN THE <Ij|Km>
C
       CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
C  (Ij,Ka) (L,R) <----- (Ij,Km) (L,L) (m,a) (X)
C
       CALL DFINDT(ICORE(ID),ICORE(IW),UIAB,DISSYD,NUMSYD,
     &             DISSYW,NUMSYW,POP(1,1),POP(1,2),VRT(1,2),
     &             IRREPL,IRREPR,IRREPX,IOFFU(1,2),2)
C
C FOR GEOMETRICAL PERTURBATIONS, INCREMENT INTEGRAL DERIVATIVE LISTS
C
       IF(STERM) THEN
C
        ID2=IW
        CALL GETLST(ICORE(ID2),1,NUMSYD,1,IRREPR,310)
        CALL SAXPY(NUMSYD*DISSYD,ONE,ICORE(ID2),1,ICORE(ID),1)
C
       ENDIF
C
       CALL PUTLST(ICORE(ID),1,NUMSYD,1,IRREPR,310)
       call checksum('D2IJKA d',icore(Id),numsyd*dissyd)
C
1000  CONTINUE
C
C ALL DONE SO FAR, WRITE INFO AND GET CPU TIMINGS
C
      CALL TIMER(1)
      write(6,6000) TIMENEW
6000  FORMAT(' Integral derivatives d <ij||ka>/d',
     &           ' chi have been formed in ',f5.1,' seconds.')
C
      RETURN
C
      END
