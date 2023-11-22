      SUBROUTINE R2L2Y1(ICORE,MAXCOR,IUHF,LISTL2,LISTL2RS,LISTR2,
     &   LISTR2RS,LISTZ,LISTZ2)
C
C DRIVER FOR CALCULATING CONTRIBUTIONS OF R2 AND L2 VECTORS
C TO THE Y1 INTERMEDIATE.
C
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ZILCH,X,SDOT
      DIMENSION ICORE(MAXCOR),I0F(2),I0Y(2)
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ZILCH/1.0D0,0.0D0/
      I0F(1)=1
      CALL GETLST(ICORE(I0F(1)),1,1,1,1,93)
      IF(IUHF.NE.0)THEN
       I0F(2)=I0F(1)+IINTFP*NT(1)
       CALL GETLST(ICORE(I0F(2)),1,1,1,2,93)
      ELSE
       I0F(2)=I0F(1)
      ENDIF
      I0Y(1)=I0F(2)+IINTFP*NT(2)
      IF(IUHF.NE.0)THEN
       I0Y(2)=I0Y(1)+IINTFP*NT(1)
      ELSE
       I0Y(2)=I0Y(1)
      ENDIF
      I000=I0Y(2)+IINTFP*NT(2)
      MXCOR=MAXCOR-I000+1
C
C INITIALIZE INTERMEDIATE WITH F(AI) * L(EF,IJ) * R(EF,IJ).
C
C L and R fully contracted, need R connected to F
      if (.false.) then
      X=ZILCH
      DO 10 ISPIN=3,3-2*IUHF,-1
       LISTL=LISTL2-1+ISPIN
       LISTR=LISTR2-1+ISPIN
       DO 20 IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        NUMDIS=IRPDPD(IRREPR,ISYTYP(2,LISTR))
        DISSIZ=IRPDPD(IRREPL,ISYTYP(1,LISTR))
        MAXT=MAX(NUMDIS,DISSIZ)
        I010=I000+IINTFP*NUMDIS*DISSIZ
        I020=I010+IINTFP*NUMDIS*DISSIZ
        CALL GETLST (ICORE(I000),1,NUMDIS,1,IRREPR,LISTL)
        CALL GETLST (ICORE(I010),1,NUMDIS,1,IRREPR,LISTR)
        IF(IUHF.EQ.0)THEN
         ITMP1=I020
         ITMP2=ITMP1+IINTFP*MAXT
         ITMP3=ITMP2+IINTFP*MAXT
         CALL SPINAD1(IRREPR,POP(1,1),DISSIZ,ICORE(I000),
     &                ICORE(ITMP1),ICORE(ITMP2))
        ENDIF
        X=X+SDOT(NUMDIS*DISSIZ,ICORE(I000),1,ICORE(I010),1)
20     CONTINUE
10    CONTINUE
      DO 11 ISPIN=1,1+IUHF
       CALL GETLST(ICORE(I0Y(ISPIN)),1,1,1,ISPIN,93)
       CALL SSCAL (NT(ISPIN),X,ICORE(I0Y(ISPIN)),1)
11    CONTINUE
      else
       DO 12 ISPIN=1,1+IUHF
        CALL ZERO(ICORE(I0Y(ISPIN)),NT(ISPIN))
12     CONTINUE
      endif
C
C PUT APPROPRIATE G INTERMEDIATES [L*R] ON LISTS 191 & 192
C
      CALL GFORMG(IRREPX,IRREPX,LISTL2,LISTR2,100,ICORE(I000),MXCOR,0,
     &   ONE,IUHF)
C
C CALCULATE CONTRACTIONS
C
C R and F fully contracted, need F to have an open line
c      CALL R2L2Y1A(ICORE(I0Y(1)),ICORE(I0F(1)),ICORE(I000),MXCOR,IUHF,
c     &   LISTL2RS,LISTR2RS)
      CALL R2L2Y1B(ICORE(I0Y(1)),ICORE(I0F(1)),ICORE(I000),MXCOR,IUHF)
      CALL R2L2Y1C(ICORE(I0Y(1)),ICORE(I0F(1)),ICORE(I000),MXCOR,IUHF)
      CALL R2L2Y1D(ICORE(I0Y(1)),ICORE(I000),MXCOR,IUHF)
      CALL R2L2Y1E(ICORE(I0Y(1)),ICORE(I000),MXCOR,IUHF)
C No open line on W
c      CALL R2L2Y1F(ICORE(I0Y(1)),ICORE(I000),MXCOR,IUHF,LISTL2RS,LISTR2)
c      CALL R2L2Y1G(ICORE(I0Y(1)),ICORE(I000),MXCOR,IUHF,LISTL2RS,LISTR2)
      IF(IUHF.EQ.0)THEN
       CALL R2L2Y1HI(ICORE(I0Y(1)),ICORE(I000),MXCOR,IUHF,LISTZ,LISTL2RS
     &      ,LISTR2RS)
      ELSE
       CALL RNGINT(ICORE(I000),MXCOR,IUHF,IRREPX,LISTL2,LISTL2RS,
     &      LISTR2,LISTR2RS,LISTZ,LISTZ2)
       CALL RNGINTRS(ICORE(I000),MXCOR,IUHF,IRREPX,1,LISTZ,LISTZ2)
       CALL R2L2Y1H(ICORE(I0Y(1)),ICORE(I000),MXCOR,IUHF,LISTZ,LISTZ2)
       CALL R2L2Y1I(ICORE(I0Y(1)),ICORE(I000),MXCOR,IUHF,LISTZ,LISTZ2)
      ENDIF
      CALL R2L2Y1J(ICORE(I0Y(1)),ICORE(I000),MXCOR,IUHF,LISTL2,LISTR2)
      CALL R2L2Y1K(ICORE(I0Y(1)),ICORE(I000),MXCOR,IUHF,LISTL2,LISTR2)
C
      DO 13 ISPIN=1,1+IUHF
       CALL GETLST(ICORE(I000),1,1,1,2+ISPIN,90)
       CALL SAXPY (NT(ISPIN),ONE,ICORE(I000),1,ICORE(I0Y(ISPIN)),1)
       CALL PUTLST(ICORE(I0Y(ISPIN)),1,1,1,2+ISPIN,90)
13    CONTINUE
C
      CALL RNGINTRS(ICORE(I000),MXCOR,IUHF,IRREPX,0,LISTZ,LISTZ2)
C
      RETURN
      END