      SUBROUTINE GETDIAG2(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE FORMS THE DIAGONAL PART OF THE TWO-BODY
C CONTRIBUTION TO THE DOUBLE-DOUBLE BLOCK OF THE EFFECTIVE
C HAMILTONIAN.  THESE QUANTITIES ARE THEN WRITTEN OUT OVER
C THE FOCK MATRIX ELEMENTS ON LISTS 91-93.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      COMMON/SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
      INDXF(I,J,N)=I+(J-1)*N
      INDXT(I,J)=I+(J*(J-1)/2)
      NNM1O2(I)=((I-1)*I)/2
C
C COLLECT ALL W(Ij,Ij) ELEMENTS 
C
      LISTW=53
      LENGTH=NOCCO(1)*NOCCO(2)
      I000=1
      I010=I000+IINTFP*LENGTH
      IOFFW=I000
      DO 100 IRREP=1,NIRREP
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
       CALL SCOPY (DISSIZ,ICORE(I010),DISSIZ+1,ICORE(IOFFW),1)
       IOFFW=IOFFW+DISSIZ*IINTFP
100   CONTINUE
      CALL PUTREC(20,'JOBARC','WDIAG53 ',IINTFP*LENGTH,ICORE)
C
      IF(IUHF.NE.0)THEN
       DO 102 ISPIN=1,2
        LISTW=50+ISPIN
        LENGTH=NNM1O2(NOCCO(ISPIN))
        I000=1
        I010=I000+IINTFP*LENGTH
        IOFFW=I000
        DO 103 IRREP=1,NIRREP
         NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
         DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
         CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
         CALL SCOPY (DISSIZ,ICORE(I010),DISSIZ+1,ICORE(IOFFW),1)
         IOFFW=IOFFW+DISSIZ*IINTFP
103     CONTINUE
        IF(ISPIN.EQ.1)THEN
         CALL PUTREC(20,'JOBARC','WDIAG51 ',IINTFP*LENGTH,ICORE)
        ELSEIF(ISPIN.EQ.2)THEN
         CALL PUTREC(20,'JOBARC','WDIAG52 ',IINTFP*LENGTH,ICORE)
        ENDIF
102    CONTINUE
      ENDIF

C
C COLLECT ALL W(Ab,Ab) ELEMENTS FOR A,b;C,d STORAGE MODE
C
      IF(ISYTYP(1,233).NE.5.AND.IFLAGS(93).NE.2)THEN
       DO 99 ISPIN=3,3-2*IUHF,-1
        IF(ISPIN.EQ.3)THEN
         LENGTH=NVRTO(1)*NVRTO(2)
        ELSEIF(ISPIN.EQ.2)THEN
         LENGTH=NNM1O2(NVRTO(ISPIN))
        ELSEIF(ISPIN.EQ.1)THEN
         LENGTH=NNM1O2(NVRTO(ISPIN))
        ENDIF
        LISTW=230+ISPIN
        I000=1
        I010=I000+IINTFP*LENGTH
        IOFFW=I000
        DO 200 IRREP=1,NIRREP
         NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
         DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
         I020=I010+IINTFP*DISSIZ*NUMDIS
         MXCOR=MAXCOR-I010+1
         IF(DISSIZ.NE.0)THEN
          NINCOR=MXCOR/(DISSIZ*IINTFP)
         ELSE
          NINCOR=NUMDIS
         ENDIF 
         NLEFT=NUMDIS
         NFIRST=1
         IOFFW0=IOFFW
1        NREAD=MIN(NLEFT,NINCOR)
         CALL GETLST(ICORE(I010),NFIRST,NREAD,1,IRREP,LISTW)
         ISTART=I010+(NFIRST-1)*IINTFP
         CALL SCOPY (NREAD,ICORE(ISTART),DISSIZ+1,ICORE(IOFFW0),1)
         IOFFW0=IOFFW0+NREAD*IINTFP
         NLEFT=NLEFT-NREAD
         NFIRST=NFIRST+NREAD
         IF(NLEFT.NE.0)GOTO 1
         IOFFW=IOFFW+DISSIZ*IINTFP
200     CONTINUE
        IF(ISPIN.EQ.3)THEN
         CALL PUTREC(20,'JOBARC','WDIAG33 ',IINTFP*LENGTH,ICORE)
        ELSEIF(ISPIN.EQ.2)THEN
         CALL PUTREC(20,'JOBARC','WDIAG32 ',IINTFP*LENGTH,ICORE)
        ELSEIF(ISPIN.EQ.1)THEN
         CALL PUTREC(20,'JOBARC','WDIAG31 ',IINTFP*LENGTH,ICORE)
        ENDIF
99     CONTINUE
C
C COLLECT ALL W(Ab,Ab) ELEMENTS FOR A<=b;C,d STORAGE MODE
C
      ELSEIF(IFLAGS(93).NE.2)THEN 
C
       LISTW=233
       LENGTH=NVRTO(1)*NVRTO(2)
       I000=1
       I010=I000+IINTFP*LENGTH
       ISTART=0
       CALL ZERO(ICORE,LENGTH)
       DO 202 IRREP=1,NIRREP
        IGET=0
        ISTART0=ISTART
        DO 203 IRREPJ=1,NIRREP
         IRREPI=DIRPRD(IRREPJ,IRREP)
         NUMJ=VRT(IRREPJ,1)
         NUMI=VRT(IRREPI,1)
         DO 204 J=1,NUMJ
          DO 205 I=1,NUMI
           IGET=IGET+1
           IF(IRREPI.LT.IRREPJ.OR.(IRREPI.EQ.IRREPJ.AND.I.LE.J))THEN
            CALL GETLST(ICORE(I010),IGET,1,1,IRREP,LISTW)
            IF(IRREPI.EQ.IRREPJ)THEN
             IMIN=MIN(I,J)
             IMAX=MAX(I,J)
             IPOS1=I010+(INDXT(IMIN,IMAX)+
     &                   ISYMOFF(IRREPJ,IRREP,5)-2)*IINTFP
            ELSE
             IPOS1=I010+(INDXF(I,J,NUMI)+
     &                   ISYMOFF(IRREPJ,IRREP,5)-2)*IINTFP
            ENDIF
            ILOC1=ISTART0+1+
     &           (INDXF(I,J,NUMI)+ISYMOFF(IRREPJ,IRREP,13)-2)*IINTFP
            ILOC2=ISTART0+1+
     &           (INDXF(J,I,NUMJ)+ISYMOFF(IRREPI,IRREP,13)-2)*IINTFP
            CALL SCOPY (1,ICORE(IPOS1),1,ICORE(ILOC1),1)
            CALL SCOPY (1,ICORE(IPOS1),1,ICORE(ILOC2),1)
           ENDIF
205       CONTINUE
204      CONTINUE
         ISTART=ISTART+NUMI*NUMJ*IINTFP
203     CONTINUE
202    CONTINUE 
       CALL PUTREC(20,'JOBARC','WDIAG33 ',IINTFP*LENGTH,ICORE)
      ELSE
       DO 104 ISPIN=3,3-2*IUHF,-1
        IF(ISPIN.EQ.3)THEN
         LENGTH=NVRTO(1)*NVRTO(2)
        ELSEIF(ISPIN.EQ.2)THEN
         LENGTH=NNM1O2(NVRTO(ISPIN))
        ELSEIF(ISPIN.EQ.1)THEN
         LENGTH=NNM1O2(NVRTO(ISPIN))
        ENDIF
        CALL ZERO(ICORE,LENGTH)
        IF(ISPIN.EQ.3)THEN
         CALL PUTREC(20,'JOBARC','WDIAG33 ',IINTFP*LENGTH,ICORE)
        ELSEIF(ISPIN.EQ.2)THEN
         CALL PUTREC(20,'JOBARC','WDIAG32 ',IINTFP*LENGTH,ICORE)
        ELSEIF(ISPIN.EQ.1)THEN
         CALL PUTREC(20,'JOBARC','WDIAG31 ',IINTFP*LENGTH,ICORE)
        ENDIF
104    CONTINUE
      ENDIF
C
C COLLECT ALL W(Ai,Ai) ELEMENTS
C
      DO 8 ISPIN=1,1+IUHF
       LISTW=57+ISPIN
       LENGTH=NVRTO(ISPIN)*NOCCO(3-ISPIN)
       I000=1
       I010=I000+IINTFP*LENGTH
       IOFFW=I000
       IOFFW0=0
       DO 300 IRREP=1,NIRREP
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
        CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
        IOFFI=I010
        DO 30 IRREPJ=1,NIRREP
         IRREPI=DIRPRD(IRREPJ,IRREP)
         NUMJ=POP(IRREPJ,3-ISPIN)
         NUMI=VRT(IRREPI,ISPIN)
         BLKSIZ=NUMI*NUMJ
         IF(BLKSIZ.NE.0)THEN
          CALL SCOPY(BLKSIZ,ICORE(IOFFI),DISSIZ+1,ICORE(IOFFW),1)
          IOFFI=IOFFI+IINTFP*(BLKSIZ*DISSIZ+BLKSIZ)
          IOFFW=IOFFW+IINTFP*BLKSIZ
          IOFFW0=IOFFW0+BLKSIZ
         ENDIF
30      CONTINUE       
300    CONTINUE
       IF(ISPIN.EQ.1)THEN
        CALL PUTREC(20,'JOBARC','WDIAG58 ',IINTFP*LENGTH,ICORE)
       ELSE
        CALL PUTREC(20,'JOBARC','WDIAG59 ',IINTFP*LENGTH,ICORE)
       ENDIF
C
C COLLECT ALL W(AI,AI) ELEMENTS
C
       LISTW=53+ISPIN
       LENGTH=NVRTO(ISPIN)*NOCCO(ISPIN)
       I000=1
       I010=I000+IINTFP*LENGTH
       IOFFW=I000
       IOFFW0=0
       DO 400 IRREP=1,NIRREP
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
        CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
        IOFFI=I010
        DO 40 IRREPJ=1,NIRREP
         IRREPI=DIRPRD(IRREPJ,IRREP)
         NUMJ=POP(IRREPJ,ISPIN)
         NUMI=VRT(IRREPI,ISPIN)
         BLKSIZ=NUMI*NUMJ
         IF(BLKSIZ.NE.0)THEN
          CALL SCOPY(BLKSIZ,ICORE(IOFFI),DISSIZ+1,ICORE(IOFFW),1)
          IOFFI=IOFFI+IINTFP*(BLKSIZ*DISSIZ+BLKSIZ)
          IOFFW=IOFFW+IINTFP*BLKSIZ
          IOFFW0=IOFFW0+BLKSIZ
         ENDIF
40      CONTINUE       
400    CONTINUE
       IF(ISPIN.EQ.1)THEN
        CALL PUTREC(20,'JOBARC','WDIAG54 ',IINTFP*LENGTH,ICORE)
       ELSE
        CALL PUTREC(20,'JOBARC','WDIAG55 ',IINTFP*LENGTH,ICORE)
       ENDIF
C
8     CONTINUE
C
      RETURN
      END