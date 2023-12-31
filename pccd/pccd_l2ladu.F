      SUBROUTINE PCCD_L2LADU(ICORE,MAXCOR,IUHF,ITYPE)
C
C DRIVER FOR LADDER CONTRIBUTIONS TO MBPT(3) ENERGY, USING SYMMETRY.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL MBPT4,CC
      INTEGER DIRPRD,DISSYW,DISSYT,POP,VRT
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
C
      MXCOR=MAXCOR
C
      MBPT4=.FALSE.
      CC=.TRUE.
C
      IF(IUHF.EQ.0) THEN
       IBOT=3
      ELSE
       IBOT=1
      ENDIF
C
      DO 1000 ISPIN=IBOT,3
C
      IF(ISPIN.LT.3)THEN
C
C    AA OR BB CASE.
C
       LISTT=ISPIN+143
       LISTZ=ISPIN+60
       IF(ITYPE.EQ.1)THEN
        LISTW=50+ISPIN
        IF(MBPT4) LISTW=10+ISPIN
       ELSE
        LISTW=230+ISPIN
       ENDIF 
       DO 100 IRREP=1,NIRREP
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        DISSYX=IRPDPD(IRREP,18+ISPIN)
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        NUMSYX=IRPDPD(IRREP,20+ISPIN)
        I000=1
        I010=I000+DISSYX*NUMSYX*IINTFP
        I020=I010+DISSYX*NUMSYX*IINTFP
        ITOP=I020+MAX(3*NT(ISPIN),DISSYT*NUMSYX)*IINTFP
        IF(ITOP.GT.MAXCOR)THEN
         CALL INSMEM('LLADAA1',I030,MAXCOR)
        ENDIF
        MAXSIZE=(MXCOR-I020)/IINTFP 
C
C CAN WE DO IT WITH THE CURRENT ALGORITHM AVAILABLE ?
C
        IF(MIN(NUMSYT,NUMSYW,DISSYT,DISSYW).NE.0)THEN
         IF(MAXSIZE.GT.DISSYW)THEN
          CALL LLADAA1(ICORE(I020),ICORE(I000),ICORE(I010),MAXSIZE,
     &                 DISSYW,NUMSYW,DISSYT,NUMSYT,LISTW,LISTT,
     &                 ITYPE,IRREP,LISTZ,ISPIN,CC,CCSD)
         ELSE
          CALL INSMEM('L2LAD',I020,MXCOR)
         ENDIF
        ENDIF
100    CONTINUE
      ELSE
C
C AB SPIN CASE.  AGAIN, FIRST TRY FOR IN-CORE ALGORITHM.
C
       LISTT=146
       IF(ITYPE.EQ.1)THEN
        LISTW=53
        IF(MBPT4) LISTW=13
       ELSE
        LISTW=233
       ENDIF
       DO 110 IRREP=1,NIRREP
C
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        I000=1
        I010=I000+MAX(NUMSYT*DISSYT,
     &                3*MAX(NUMSYT,DISSYT,NT(1)))*IINTFP
        I020=I010+NUMSYT*DISSYT*IINTFP
        MAXSIZE=(MXCOR-I020)/IINTFP
C
C CAN WE DO IT WITH THE MEMORY AVAILABLE
C
        IF(MIN(NUMSYT,NUMSYW,DISSYT,DISSYW).NE.0)THEN
         IF(MAXSIZE.GT.DISSYW)THEN
          IF(IUHF.NE.0.OR.ITYPE.EQ.1)THEN
           CALL LLADAB1(ICORE(I020),ICORE(I000),ICORE(I010),MAXSIZE,
     &                  DISSYW,NUMSYW,DISSYT,NUMSYT,LISTW,LISTT,ITYPE,
     &                  IRREP,CC,CCSD)
          ELSE
           CALL LLADAB2(ICORE(I020),ICORE(I000),ICORE(I010),MAXSIZE,
     &                  DISSYW,NUMSYW,DISSYT,NUMSYT,LISTW,LISTT,ITYPE,
     &                  IRREP,CC,CCSD)
          ENDIF
C
C NO OTHER ALGORITHM AVAILABLE AND PLANNED !
C
         ELSE
          CALL INSMEM('L2LAD',I020,MXCOR)
         ENDIF
        ENDIF
110    CONTINUE
      ENDIF
1000  CONTINUE
      RETURN
      END
