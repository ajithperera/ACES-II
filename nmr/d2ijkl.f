      SUBROUTINE D2IJKL(UAIA,UAIB,SIJA,SIJB,ICORE,MXCOR,IRREPX,
     &                  IPERT,IUHF,STERM,ANTI,GRAD1)
C
C  THIS ROUTINE CALCULATES THE DERIVATIVE INTEGRALS d <Ij||Kl> / d chi
C
C   THERE ARE FOUR DIFFERENT TERMS
C
CEND
C
C CODED AUG/91  JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYW,DISSYD,POP,VRT
      LOGICAL FIELD,GEOM,ROHF,QRHF,SEMI,STERM,ANTI
      DIMENSION ICORE(MXCOR) 
      DIMENSION UAIA(1),UAIB(1),SIJA(1),SIJB(1)
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
C  INCORE ALGORITHM
C
C  DO HERE ALPHA,BETA SPIN CASE
C
C  LOOP OVER ALL IRREPS OF THE DERIVATIVE INTEGRALS
C
      DO 1000 IRREPR=1,NIRREP
C
C  THE IRREP ON THE LEFT SIDE IS THEN GIVEN AS THE DIRECT PRODUCT
C   OF IRREPP AND IRREPR
C
       IRREPL=DIRPRD(IRREPR,IRREPX)
C
C DETERMINE LENGTH OF TARGET ARRAY
C
       NUMSYD=IRPDPD(IRREPR,ISYTYP(2,13))
       DISSYD=IRPDPD(IRREPL,ISYTYP(1,13))
C
C ALLOCATE CORE FOR TARGET ARRAY
C
       ID=1
       IREST=ID+IINTFP*NUMSYD*DISSYD
C
C  ZERO TARGET LIST
C
        CALL ZERO(ICORE(ID),NUMSYD*DISSYD)
C
C CONSIDER FIRST ALL TERMS WHICH INCLUDE U(AI)
C
C    <Ij||Ke> U(el) + <Ij||El> U(EK)
C
C  THE IRREP ON THE LEFT SIDE IS UNCHANGED, SO IRREPL DETERMINES
C  THE LISTS OF ORIGINAL ZEROTH ORDER INTEGRALS
C
       LISTW=10
       NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
       DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
C  ALLOCATE CORE MEMORY
C
       IW=IREST
       IEND=IW+IINTFP*NUMSYW*DISSYW
C
       IF(IEND.GE.MXCOR) CALL ERREX
C
C  READ IN THE <Ij||Ke> 
C
        CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
C  (Ij;Kl)(L,R) <---- (Ij;Ke)(L,L) (el) (X)
C
C 03/08/97 missing comma was added AP
C
        CALL DFINDT(ICORE(ID),ICORE(IW),UAIB,DISSYD,NUMSYD,
     &              DISSYW,NUMSYW,POP(1,1),VRT(1,2),POP(1,2),
     &              IRREPL,IRREPR,IRREPX,IOFFU(1,2),1)
C
C  NOW DO <Ij||El> U(E,K) , BUT ONLY FOR UHF
C
        IF(IUHF.EQ.1) THEN
C
C  TRANSPOSE THE TARGET LIST    Ij;Kl --> Ij;lK
C
        ITMP=IREST
        IEND=ITMP+3*IINTFP*MAX(NUMSYD,DISSYD)
        IF(IEND.GE.MXCOR) CALL ERREX
C
        CALL SYMTR1(IRREPR,POP(1,1),POP(1,2),DISSYD,ICORE(ID),
     &              ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYD),
     &              ICORE(ITMP+2*IINTFP*DISSYD))
C
         LISTW=9
         NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
C  ALLOCATE CORE MEMORY
C
         IW=IREST
         ITMP=IW+IINTFP*NUMSYW*DISSYW
         IEND=ITMP+3*IINTFP*MAX(NUMSYW,DISSYW,DISSYD,NUMSYD)
C
         IF(IEND.GE.MXCOR) CALL ERREX
C
C  READ IN THE <Ij||El>    
C
         CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
C TRANSPOSE  Ij;El ---> Ij;lE 
C
         CALL SYMTR1(IRREPL,VRT(1,1),POP(1,2),DISSYW,ICORE(IW),
     &               ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &               ICORE(ITMP+2*IINTFP*DISSYD))
C
C  (Ij;lK) (L,R) <----- (Ij;lE) (L,L) (E,K) (X)
C
         CALL DFINDT(ICORE(ID),ICORE(IW),UAIA,DISSYD,NUMSYD,
     &               DISSYW,NUMSYW,POP(1,2),VRT(1,1),POP(1,1),
     &               IRREPL,IRREPR,IRREPX,IOFFU(1,1),1)
        ENDIF 
C
        IF(STERM) THEN
C
C CONSIDER  TERMS WHICH INCLUDE S(AE)
C
C    <Ij||Km> (-1/2 S(ml)) + <Ij||Ml> (-1/2 S(MK))
C
C  THE IRREP ON THE LEFT SIDE IS UNCHANGED, SO IRREPL DETERMINES
C  THE LISTS OF ORIGINAL ZEROTH ORDER INTEGRALS
C
         LISTW=13
         NUMSYW=IRPDPD(IRREPL,ISYTYP(2,LISTW))
         DISSYW=IRPDPD(IRREPL,ISYTYP(1,LISTW))
C
C  ALLOCATE CORE MEMORY
C
         IW=IREST
         ITMP=IW+IINTFP*NUMSYW*DISSYW
         IEND=ITMP+3*IINTFP*MAX(NUMSYW,NUMSYD,DISSYW,DISSYD)*IUHF
C
         IF(IEND.GE.MXCOR) CALL ERREX
C
C  READ IN THE <Ij||Kl> INTEGRALS
C
         CALL GETLST(ICORE(IW),1,NUMSYW,1,IRREPL,LISTW)
C
         IF(IUHF.EQ.1) THEN
C
C  TRANSPOSE THE LAST TWO INDICES :   Ij,Ml --> Ij;lM 
C
         CALL SYMTR1(IRREPL,POP(1,1),POP(1,2),DISSYW,ICORE(IW),
     &               ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &               ICORE(ITMP+2*IINTFP*DISSYW))
C
C  (Ij;lK) (L,R) <----- (Ij;lM) (L,L) (M,K) (X)
C
         CALL DFINDT(ICORE(ID),ICORE(IW),SIJA,DISSYD,NUMSYD,
     &               DISSYW,NUMSYW,POP(1,2),POP(1,1),POP(1,1),
     &               IRREPL,IRREPR,IRREPX,IOFFS1(1,1),1)
C
C  TRANSPOSE THE TARGET LIST   Ij;lK --> Ij;Kl
C
         CALL SYMTR1(IRREPR,POP(1,2),POP(1,1),DISSYD,ICORE(ID),
     &               ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYD),
     &               ICORE(ITMP+2*IINTFP*DISSYD))
C
C  TRANSPOSE THE INTEGRAL LIST
C
          CALL SYMTR1(IRREPL,POP(1,2),POP(1,1),DISSYW,ICORE(IW),
     &                ICORE(ITMP),ICORE(ITMP+IINTFP*DISSYW),
     &                ICORE(ITMP+2*IINTFP*DISSYW))
C
          ENDIF
C
C   (Ij,Kl) (L,R) <----- (Ij,Km) (L,L) (m,l) (X)
C
          CALL DFINDT(ICORE(ID),ICORE(IW),SIJB,DISSYD,NUMSYD,
     &                DISSYW,NUMSYW,POP(1,1),POP(1,2),POP(1,2),
     &                IRREPL,IRREPR,IRREPX,IOFFS1(1,2),1)
C
        ENDIF 
C
C  SYMMETRIZE THE INTEGRAL DERIVATIVES
C
c       CALL DSYMMET1(IRREPX,ICORE(ID),IRPDPD(1,ISYTYP(1,13)))
C
C  FOR RHF, DO NOW ALL TRANSPOSITIONS
C
       IF(IUHF.EQ.0) THEN
C
        ITMP=IREST
C
        CALL DSYMRHF(IRREPL,IRREPR,VRT(1,1),POP(1,1),DISSYD,
     &               ICORE(ID),ICORE(ITMP),
     &               ICORE(ITMP+IINTFP*DISSYD),
     &               ICORE(ITMP+2*IINTFP*DISSYD))
C
       ENDIF
C
       call checksum('D2IJKL',icore(Id),numsyd*dissyd)

c        if(sterm)then
c        call getlst(icore(iend),(ipert-1)*numsyd+1,numsyd,
c     &                2,1,313)
c        open(unit=16,form='formatted')
c        read(16,'((3f20.10))')(icore(iend-1+j),j=1,numsyd*dissyd)
c        call saxpy(numsyd*dissyd,one,icore(iend),1,icore(id),1)
c        call scopy(numsyd*dissyd,icore(iend),1,icore(id),1)
c        close(unit=16,status='keep')
c        endif

c       CALL PUTLST(ICORE(ID),(IPERT-1)*NUMSYD+1,NUMSYD,1,
c     &             IRREPR,313)
1000  CONTINUE
C
C ALL DONE SO FAR
C
      RETURN
      END
