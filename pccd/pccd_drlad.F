      SUBROUTINE PCCD_DRLAD(ICORE,MAXCOR,IUHF,ITYPE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYW,DISSYT
      INTEGER POP,VRT
      CHARACTER*4 MODIF
      CHARACTER*2 SPCASE(3)

      DIMENSION ICORE(MAXCOR)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON /FLAGS/ IFLAGS(100)
C
      LISTT=46
      IF(ITYPE.EQ.1)THEN
         LISTW=53
      ELSE
         LISTW=233
      ENDIF

      DO 110 IRREP=1,NIRREP
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))

        I000=1
        I010=I000+NUMSYT*DISSYT*IINTFP
        I020=I010+NUMSYT*DISSYT*IINTFP  
        MAXSIZE=(MAXCOR-I020)/IINTFP

        IF(MIN(NUMSYT,NUMSYW,DISSYT,DISSYW).NE.0)THEN
         IF(MAXSIZE.GT.DISSYW) THEN
          CALL PCCD_LADAB1(ICORE(I020),ICORE(I000),ICORE(I010),MAXSIZE,
     &                     DISSYW,NUMSYW,DISSYT,NUMSYT,LISTW,LISTT,
     &                     POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),ITYPE,
     &                     IRREP,IUHF)
         ELSE
          CALL INSMEM('PCCD_DRLAD',DISSYW,MAXSIZE)
         ENDIF
        ENDIF
110    CONTINUE

      RETURN
      END
