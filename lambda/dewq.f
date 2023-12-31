      SUBROUTINE DEWQ(ICORE,MAXCOR,IUHF)
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
C COMPUTE AB CONTRIBUTION TO E4D.
C
      DO 10 IRREP=1,NIRREP
       NDSSYM=IRPDPD(IRREP,ISYTYP(2,63))
       DISSYM=IRPDPD(IRREP,ISYTYP(1,63))
       NSIZE=DISSYM*NDSSYM
       I000=1
       I010=I000+IINTFP*NSIZE
       I020=I010+IINTFP*NSIZE
       CALL GETLST(ICORE(I000),1,NDSSYM,1,IRREP,63)
       CALL GETLST(ICORE(I010),1,NDSSYM,1,IRREP,50)
       CALL DEWQ1(ICORE(I000),ICORE(I010),NSIZE)
       CALL PUTLST(ICORE(I000),1,NDSSYM,1,IRREP,63)
10    CONTINUE
C
C COMPUTE AA AND BB (BB FOR UHF ONLY) CONTRIBUTIONS TO E4D.
C
      DO 100 ISPIN=1,1+IUHF
       NLIST=60+ISPIN
       DO 20 IRREP=1,NIRREP
        NDSSYM=IRPDPD(IRREP,ISYTYP(2,NLIST))
        DISSYM=IRPDPD(IRREP,ISYTYP(1,NLIST))
        NSIZE=DISSYM*NDSSYM
        I000=1
        I010=I000+IINTFP*NSIZE
        I020=I010+IINTFP*NSIZE
        CALL GETLST(ICORE(I000),1,NDSSYM,1,IRREP,NLIST)
        CALL GETLST(ICORE(I010),1,NDSSYM,1,IRREP,47+ISPIN)
        CALL DEWQ1(ICORE(I000),ICORE(I010),NSIZE)
        CALL PUTLST(ICORE(I000),1,NDSSYM,1,IRREP,NLIST)
20     CONTINUE
100   CONTINUE
      RETURN
      END
