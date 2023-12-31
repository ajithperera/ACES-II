      SUBROUTINE R2L1Y1(ICORE,MAXCOR,IUHF)
C
C DRIVER FOR CALCULATING CONTRIBUTIONS OF R2 AND L1 VECTORS
C TO THE Y1 INTERMEDIATE.
C  
CEND
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ZILCH,X,SDOT
      DIMENSION ICORE(MAXCOR),I0Y(2),I0L(2)
      LOGICAL MBPT2,CC,CCD,RPA,DRPA,LCCD,LCCSD,CC2
      COMMON/REFTYPE/MBPT2,CC,CCD,RPA,DRPA,LCCD,LCCSD,CC2
      COMMON/STATSYM/IRREPX
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA ONE,ZILCH/1.0D0,0.0D0/
C
      I0L(1)=1
      CALL GETLST(ICORE(I0L(1)),1,1,1,1,490)
      IF(IUHF.NE.0)THEN
       I0L(2)=I0L(1)+IINTFP*IRPDPD(IRREPX,9)
       CALL GETLST(ICORE(I0L(2)),1,1,1,2,490)
      ELSE
       I0L(2)=I0L(1)
      ENDIF 
      I0Y(1)=I0L(2)+IINTFP*IRPDPD(IRREPX,10)
      IF(IUHF.NE.0)THEN
       I0Y(2)=I0Y(1)+IINTFP*NT(1)
      ELSE
       I0Y(2)=I0Y(1)
      ENDIF
      I000=I0Y(2)+IINTFP*NT(2)
      MXCOR=MAXCOR-I000+1
C
C INITIALIZE INTERMEDIATE WITH L(AI) * W(EF,IJ) * R(EF,IJ) FOR
C TOTALLY SYMMETRIC EOM VECTORS ONLY
C
      IF(IRREPX.EQ.1 .AND. .NOT. (MBPT2 .OR. CC2))THEN
       X=ZILCH 
       DO 10 ISPIN=3,3-2*IUHF,-1
        LISTW=13+ISPIN
        LISTR=460+ISPIN
        DO 20 IRREP=1,NIRREP
         NUMDIS=IRPDPD(IRREP,ISYTYP(2,LISTW))
         DISSIZ=IRPDPD(IRREP,ISYTYP(1,LISTW))
         MAXT=MAX(NUMDIS,DISSIZ)
         I010=I000+IINTFP*NUMDIS*DISSIZ
         I020=I010+IINTFP*NUMDIS*DISSIZ
         CALL GETLST (ICORE(I000),1,NUMDIS,1,IRREP,LISTW)
         CALL GETLST (ICORE(I010),1,NUMDIS,1,IRREP,LISTR)
         IF(IUHF.EQ.0)THEN
          ITMP1=I020
          ITMP2=ITMP1+IINTFP*MAXT  
          ITMP3=ITMP2+IINTFP*MAXT
          IEND=ITMP3+IINTFP*MAXT
          IF(IEND.GE.MAXCOR) CALL INSMEM('R2L1Y1',IEND,MAXCOR)
          CALL SPINAD1(IRREP,POP(1,1),DISSIZ,ICORE(I000),
     &                 ICORE(ITMP1),ICORE(ITMP2))
         ENDIF
         X=X+SDOT(NUMDIS*DISSIZ,ICORE(I000),1,ICORE(I010),1)
20      CONTINUE
10     CONTINUE
       DO 11 ISPIN=1,1+IUHF
        CALL GETLST(ICORE(I0Y(ISPIN)),1,1,1,ISPIN,490)
        CALL SSCAL (NT(ISPIN),X,ICORE(I0Y(ISPIN)),1)
11     CONTINUE
      ELSE
       DO 12 ISPIN=1,1+IUHF
        CALL ZERO(ICORE(I0Y(ISPIN)),NT(ISPIN))
12     CONTINUE
      ENDIF 
C
C                                    
C PUT APPROPRIATE G INTERMEDIATES [W*R] ON LISTS 491 & 492
C
      CALL GFORMG(IRREPX,1,461,14,400,ICORE(I000),MXCOR,0,ONE,ONE,
     &            IUHF)
C
C CALCULATE CONTRACTIONS
C
      CALL R2L1Y1A(ICORE(I0Y(1)),ICORE(I0L(1)),ICORE(I000),MXCOR,IUHF)
      CALL R2L1Y1B(ICORE(I0Y(1)),ICORE(I0L(1)),ICORE(I000),MXCOR,IUHF)
      CALL R2L1Y1C(ICORE(I0Y(1)),ICORE(I0L(1)),ICORE(I000),MXCOR,IUHF)
C
      DO 13 ISPIN=1,1+IUHF
       CALL PUTLST(ICORE(I0Y(ISPIN)),1,1,1,2+ISPIN,90)
13    CONTINUE
C
      RETURN
      END
