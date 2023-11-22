      SUBROUTINE GETDIAG(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE FORMS THE DIAGONAL PART OF THE TWO-BODY
C CONTRIBUTION TO THE DOUBLE-DOUBLE BLOCK OF THE EFFECTIVE
C HAMILTONIAN.  THESE QUANTITIES ARE THEN WRITTEN OUT OVER
C THE FOCK MATRIX ELEMENTS ON LISTS 91-93.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL CIS,EOMCC,FULDIAG,INCORE,READGUES
      DOUBLE PRECISION ONE,ONEM,ZILCH,FACT,SNRM2
      CHARACTER*8 LABEL(2),LABEL2(2)
      DIMENSION ICORE(MAXCOR),IBOT(8)
      COMMON/SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/METH/CIS,EOMCC,FULDIAG,INCORE,READGUES
      COMMON/DIAGW/IOFFIJ(8,8),IOFFAB(8,8),IOFFAI(8,8)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
C
C COLLECT ALL W(IJ,IJ) ELEMENTS 
C
      LISTW=53
      LENGTH=NOCCO(1)*NOCCO(2)
      I000=1
      I010=I000+IINTFP*LENGTH
      IOFFW=I000
      IOFFIJ0=0
      DO 100 IRREP=1,NIRREP
       NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
       DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
       CALL GETLST(ICORE(I010),1,NUMDIS,1,IRREP,LISTW)
       CALL SCOPY (DISSIZ,ICORE(I010),DISSIZ+1,ICORE(IOFFW),1)
       DO 101 IRREPJ=1,NIRREP
        IOFFIJ(IRREPJ,IRREP)=IOFFIJ0+ISYMOFF(IRREPJ,IRREP,14)-1
101    CONTINUE
       IOFFW=IOFFW+DISSIZ*IINTFP
       IOFFIJ0=IOFFIJ0+DISSIZ
100   CONTINUE
      CALL PUTREC(20,'JOBARC','DIAGWIJ3',IINTFP*LENGTH,ICORE)
C
C COLLECT ALL W(AB,AB) ELEMENTS
C
      LISTW=233
      LENGTH=NVRTO(1)*NVRTO(2)
      I000=1
      I010=I000+IINTFP*LENGTH
      IOFFW=I000
      IOFFAB0=0
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
1      NREAD=MIN(NLEFT,NINCOR)
       CALL GETLST(ICORE(I010),NFIRST,NREAD,1,IRREP,LISTW)
       ISTART=I010+(NFIRST-1)*IINTFP
       CALL SCOPY (NREAD,ICORE(ISTART),DISSIZ+1,ICORE(IOFFW0),1)
       IOFFW0=IOFFW0+NREAD*IINTFP
       NLEFT=NLEFT-NREAD
       NFIRST=NFIRST+NREAD
       IF(NLEFT.NE.0)GOTO 1
       DO 201 IRREPJ=1,NIRREP
        IOFFAB(IRREPJ,IRREP)=IOFFAB0+ISYMOFF(IRREPJ,IRREP,13)-1
201    CONTINUE
       IOFFAB0=IOFFAB0+DISSIZ
       IOFFW=IOFFW+DISSIZ*IINTFP
200   CONTINUE
      CALL PUTREC(20,'JOBARC','DIAGWAB3',IINTFP*LENGTH,ICORE)
C
C COLLECT ALL W(AI,AI) ELEMENTS
C
      LISTW=58
      LENGTH=NVRTO(1)*NOCCO(2)
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
        NUMJ=POP(IRREPJ,2)
        NUMI=VRT(IRREPI,1)
        BLKSIZ=NUMI*NUMJ
        IOFFAI(IRREPJ,IRREP)=IOFFW0
        IF(BLKSIZ.NE.0)THEN
         CALL SCOPY(BLKSIZ,ICORE(IOFFI),DISSIZ+1,ICORE(IOFFW),1)
         IOFFI=IOFFI+IINTFP*(BLKSIZ*DISSIZ+BLKSIZ)
         IOFFW=IOFFW+IINTFP*BLKSIZ
         IOFFW0=IOFFW0+BLKSIZ
        ENDIF
30     CONTINUE       
300   CONTINUE
      CALL PUTREC(20,'JOBARC','DIAGWRN1',IINTFP*LENGTH,ICORE)
C
C COLLECT ALL W(AI,AI) ELEMENTS
C
      LISTW=54
      LENGTH=NVRTO(1)*NOCCO(2)
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
        NUMJ=POP(IRREPJ,1)
        NUMI=VRT(IRREPI,1)
        BLKSIZ=NUMI*NUMJ
        IOFFAI(IRREPJ,IRREP)=IOFFW0
        IF(BLKSIZ.NE.0)THEN
         CALL SCOPY(BLKSIZ,ICORE(IOFFI),DISSIZ+1,ICORE(IOFFW),1)
         IOFFI=IOFFI+IINTFP*(BLKSIZ*DISSIZ+BLKSIZ)
         IOFFW=IOFFW+IINTFP*BLKSIZ
         IOFFW0=IOFFW0+BLKSIZ
        ENDIF
40     CONTINUE       
400   CONTINUE
      CALL PUTREC(20,'JOBARC','DIAGWRN2',IINTFP*LENGTH,ICORE)
C
      RETURN
      END
