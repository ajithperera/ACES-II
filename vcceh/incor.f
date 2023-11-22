C
C THIS ROUTINE SETS CERTAIN LISTS TO EXIST IN THE AUXILIARY CACHE, SO
C  THAT IT IS NEVER NECESSARY TO GO TO DISK TO GET THESE QUANTITIES.
C
      SUBROUTINE INCOR(I0,IUHF)
      IMPLICIT INTEGER (A-Z)
      COMMON / / ICORE(1)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /IOSTAT/ NOPBUF,NOPDSK
      COMMON /ISTART/ IGNORE,ICRSIZ
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /AUXCACHE/ QUIKGET(10,500)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
C
      print *,'@INCOR: vcceh has not been modified to use INCORE safely'
      print *,'        The calculation will proceed with INCORE=OFF.'
      return

      CALL IZERO(QUIKGET,5000)
C
      IF(IFLAGS(2).LE.1) RETURN
      IF(IFLAGS(35).EQ.0)RETURN
C
      IF(IFLAGS(35).EQ.1)THEN
C
C HOLD EVERYTHING BUT ABCD INTEGRALS IN CORE!
C
       DO 110 ILIST=1,199
        DO 111 ISUBLT=1,10
         NUMDIS=aces_list_cols(ISUBLT,ILIST)
         DISSIZ=aces_list_rows(ISUBLT,ILIST)
         IF(NUMDIS.NE.0)THEN
          CALL GETLST(ICORE(I0),1,NUMDIS,1,ISUBLT,ILIST)
          QUIKGET(ISUBLT,ILIST)=I0
          I0=I0+NUMDIS*DISSIZ*IINTFP
          ICRSIZ=ICRSIZ-NUMDIS*DISSIZ*IINTFP
         ENDIF
111     CONTINUE
110    CONTINUE
C
      ELSEIF(IFLAGS(35).EQ.6)THEN
C
C HOLD ABSOLUTELY EVERYTHING IN CORE
C
       DO 210 ILIST=1,500
        DO 211 ISUBLT=1,10
         NUMDIS=aces_list_cols(ISUBLT,ILIST)
         DISSIZ=aces_list_rows(ISUBLT,ILIST)
         IF(NUMDIS.NE.0)THEN
          CALL GETLST(ICORE(I0),1,NUMDIS,1,ISUBLT,ILIST)
          QUIKGET(ISUBLT,ILIST)=I0
          I0=I0+NUMDIS*DISSIZ*IINTFP
          ICRSIZ=ICRSIZ-NUMDIS*DISSIZ*IINTFP
         ENDIF
211     CONTINUE
210    CONTINUE
C
      ELSEIF(IFLAGS(35).EQ.5)THEN
C
C HOLD IJKL, IJKA AND SINGLES IN CORE
C
       DO 212 ILIST=1,199
        J=ILIST
        IF(J.GE.7.AND.J.LE.13.OR.J.GE.51.AND.J.LE.53.OR.
     &     J.GE.90.AND.J.LE.94)THEN
         DO 213 ISUBLT=1,10
          NUMDIS=aces_list_cols(ISUBLT,ILIST)
          DISSIZ=aces_list_rows(ISUBLT,ILIST)
          IF(NUMDIS.NE.0)THEN
           CALL GETLST(ICORE(I0),1,NUMDIS,1,ISUBLT,ILIST)
           QUIKGET(ISUBLT,ILIST)=I0
           I0=I0+NUMDIS*DISSIZ*IINTFP
           ICRSIZ=ICRSIZ-NUMDIS*DISSIZ*IINTFP
          ENDIF
213      CONTINUE
        ENDIF
212    CONTINUE
C
      ELSEIF(IFLAGS(35).EQ.4)THEN
C
C HOLD EVERYTHING BUT ABCD AND ABCI INTEGRALS IN CORE!
C
       DO 214 ILIST=1,199
        IF(MOD(ILIST,100).LT.27.OR.MOD(ILIST,100).GT.30)THEN
         DO 215 ISUBLT=1,10
          NUMDIS=aces_list_cols(ISUBLT,ILIST)
          DISSIZ=aces_list_rows(ISUBLT,ILIST)
          IF(NUMDIS.NE.0)THEN
           CALL GETLST(ICORE(I0),1,NUMDIS,1,ISUBLT,ILIST)
           QUIKGET(ISUBLT,ILIST)=I0
           I0=I0+NUMDIS*DISSIZ*IINTFP
           ICRSIZ=ICRSIZ-NUMDIS*DISSIZ*IINTFP
          ENDIF
215      CONTINUE
        ENDIF
214    CONTINUE
C
      ELSEIF(IFLAGS(35).EQ.2)THEN
C
C HOLD T2 AND INCREMENTS IN CORE
C
      DO 10 ILIST=444+2*(1-IUHF),446,1
       DO 11 IRREP=1,NIRREP
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,ILIST))
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,ILIST))
        CALL GETLST(ICORE(I0),1,NUMDIS,1,IRREP,ILIST)
        QUIKGET(IRREP,ILIST)=I0
        I0=I0+NUMDIS*DISSIZ*IINTFP
        ICRSIZ=ICRSIZ-NUMDIS*DISSIZ*IINTFP
11     CONTINUE
10    CONTINUE
C
      DO 20 ILIST=461+2*(1-IUHF),463,1
       DO 21 IRREP=1,NIRREP
        NUMDIS=IRPDPD(IRREP,ISYTYP(2,ILIST))
        DISSIZ=IRPDPD(IRREP,ISYTYP(1,ILIST))
        CALL GETLST(ICORE(I0),1,NUMDIS,1,IRREP,ILIST)
        QUIKGET(IRREP,ILIST)=I0
        I0=I0+NUMDIS*DISSIZ*IINTFP
        ICRSIZ=ICRSIZ-NUMDIS*DISSIZ*IINTFP
21     CONTINUE
20    CONTINUE
C
C FOR QCISD AND CCSD KEEP THE SINGLE STUFF IN CORE AS WELL
C
      IF(IFLAGS(2).GE.10) THEN
C
       DO 30 ISPIN=1,IUHF+1
        CALL GETLST(ICORE(I0),1,1,1,ISPIN,490)
        QUIKGET(ISPIN,90)=I0
        I0=I0+NT(ISPIN)*IINTFP
        ICRSIZ=ICRSIZ-IINTFP*NT(ISPIN)
        CALL GETLST(ICORE(I0),1,1,1,ISPIN+2,490)
        QUIKGET(ISPIN+2,90)=I0
        I0=I0+NT(ISPIN)*IINTFP
        ICRSIZ=ICRSIZ-IINTFP*NT(ISPIN)
        CALL GETLST(ICORE(I0),1,1,1,ISPIN,491)
        QUIKGET(ISPIN,91)=I0
        I0=I0+NFMI(ISPIN)*IINTFP
        ICRSIZ=ICRSIZ-IINTFP*NFMI(ISPIN)
        CALL GETLST(ICORE(I0),1,1,1,ISPIN,492)
        QUIKGET(ISPIN,92)=I0
        I0=I0+NFEA(ISPIN)*IINTFP
        ICRSIZ=ICRSIZ-IINTFP*NFEA(ISPIN)
        CALL GETLST(ICORE(I0),1,1,1,ISPIN,493)
        QUIKGET(ISPIN,93)=I0
        I0=I0+NT(ISPIN)*IINTFP
        ICRSIZ=ICRSIZ-IINTFP*NT(ISPIN)
30     CONTINUE 
C
       IF(IFLAGS(11).EQ.2.OR.IFLAGS(77).NE.0)THEN
        DO 40 ISPIN=1,IUHF+1
        CALL GETLST(ICORE(I0),1,1,1,ISPIN+2,491)
        QUIKGET(ISPIN+2,91)=I0
        I0=I0+NFMI(ISPIN)*IINTFP
        ICRSIZ=ICRSIZ-IINTFP*NFMI(ISPIN)
        CALL GETLST(ICORE(I0),1,1,1,ISPIN+2,492)
        QUIKGET(ISPIN+2,92)=I0
        I0=I0+NFEA(ISPIN)*IINTFP
        ICRSIZ=ICRSIZ-IINTFP*NFEA(ISPIN)
        CALL GETLST(ICORE(I0),1,1,1,ISPIN+2,493)
        QUIKGET(ISPIN+2,93)=I0
        I0=I0+NT(ISPIN)*IINTFP
        ICRSIZ=ICRSIZ-IINTFP*NT(ISPIN)
40      CONTINUE
       ENDIF
      ENDIF
      ENDIF
      RETURN
      END