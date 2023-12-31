       SUBROUTINE PCCD_QUAD1(ICORE,MAXCOR,IUHF,FACT,LCCD)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYW,DISSYT,POP1,POP2,VRT1,VRT2
      DIMENSION ICORE(MAXCOR)
      LOGICAL LCCD

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYM/ POP1(8),POP2(8),VRT1(8),VRT2(8),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)

      LISTT=46
      LISTW=16

      DO 110 IRREP=1,NIRREP
          DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
          DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
          NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
          NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
          I001=1
          I002=I001+IINTFP*NUMSYW*DISSYW
          I003=I002+IINTFP*NUMSYT*DISSYT

          IF(MIN(NUMSYT,NUMSYW,DISSYT,DISSYW).NE.0)THEN
             I004=I003+IINTFP*NUMSYT*NUMSYW
             IF(I004.LT.MAXCOR) THEN
                CALL PCCD_Q1ALL(ICORE(I001),ICORE(I002),ICORE(I003),
     &                          DISSYW,DISSYT,NUMSYW,NUMSYT,FACT,
     &                          IRREP,LCCD)
             ELSE
                CALL INSMEM("PCCD_QUAD1",I004,MAXCOR)
             ENDIF
          ELSE
          IF (NUMSYW.NE.0) THEN
             NLIST1=13
             NLIST2=53
             CALL GETLST(ICORE(I001),1,NUMSYW,2,IRREP,NLIST1)
             CALL PUTLST(ICORE(I001),1,NUMSYW,1,IRREP,NLIST2)
          ENDIF
          ENDIF 
110   CONTINUE

      RETURN
      END
