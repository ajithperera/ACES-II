      SUBROUTINE DEWQ(ICORE,MAXCOR,IUHF,DVORML)
      IMPLICIT INTEGER (A-Z)
      DIMENSION ICORE(MAXCOR)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      IF(DVORML.EQ.0)THEN
C
      WRITE(6,1000)
 1000 FORMAT(' @DEWQ-I, T ---> DT ')
C
C     ADD ENERGY DIFFERENCE D(IJ,AB) = E(I) + E(J) - E(A) - E(B) TO
C     T2 AMPLITUDE INCREMENTS (T ---> DT) ON LISTS 61-63.
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
C
      ELSE
C
      WRITE(6,1010)
 1010 FORMAT(' @DEWQ-I, DT ---> T ')
C
C     REMOVE ENERGY DIFFERENCE D(IJ,AB) = E(I) + E(J) - E(A) - E(B)
C     FROM D2T2 INCREMENTS (DT ---> T) ON LISTS 61-63.
C
      DO 110 IRREP=1,NIRREP
       NDSSYM=IRPDPD(IRREP,ISYTYP(2,63))
       DISSYM=IRPDPD(IRREP,ISYTYP(1,63))
       NSIZE=DISSYM*NDSSYM
       I000=1
       I010=I000+IINTFP*NSIZE
       I020=I010+IINTFP*NSIZE
       CALL GETLST(ICORE(I000),1,NDSSYM,1,IRREP,63)
       CALL GETLST(ICORE(I010),1,NDSSYM,1,IRREP,50)
       CALL VECPRD(ICORE(I000),ICORE(I010),ICORE(I000),NSIZE)
C      CALL DEWQ1(ICORE(I000),ICORE(I010),NSIZE)
       CALL PUTLST(ICORE(I000),1,NDSSYM,1,IRREP,63)
  110 CONTINUE
C
      DO 130 ISPIN=1,1+IUHF
       NLIST=60+ISPIN
       DO 120 IRREP=1,NIRREP
        NDSSYM=IRPDPD(IRREP,ISYTYP(2,NLIST))
        DISSYM=IRPDPD(IRREP,ISYTYP(1,NLIST))
        NSIZE=DISSYM*NDSSYM
        I000=1
        I010=I000+IINTFP*NSIZE
        I020=I010+IINTFP*NSIZE
        CALL GETLST(ICORE(I000),1,NDSSYM,1,IRREP,NLIST)
        CALL GETLST(ICORE(I010),1,NDSSYM,1,IRREP,47+ISPIN)
        CALL VECPRD(ICORE(I000),ICORE(I010),ICORE(I000),NSIZE)
C       CALL DEWQ1(ICORE(I000),ICORE(I010),NSIZE)
        CALL PUTLST(ICORE(I000),1,NDSSYM,1,IRREP,NLIST)
  120  CONTINUE
  130 CONTINUE
C
      ENDIF
      RETURN
      END
