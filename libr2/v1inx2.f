      SUBROUTINE V1INX2(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE FORMS THE FOLLWING CONTRIBUTION TO X(IJ,AB)
C
C    X(IJ,AB) =    1/2 SUM M,N    TAU(MN,AB) V(IJ,MN)
C
CEND
C
C  CODED AUGUST/90  JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TAU
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      INTEGER DIRPRD,DISSYW,DISSYT
      INTEGER POP,VRT
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /METH/ MBPT2,MBPt3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
C
      MXCOR=MAXCOR
      TAU=.FALSE.
      IF(CCSD) THEN
C
C    ALLOCATE CORE MEMORY FOR T1 AMPLITUDES
C
      I0TA=MXCOR+1-NTAA*IINTFP
      MXCOR=MXCOR-NTAA*IINTFP
      IF(IUHF.EQ.0) THEN
       I0TB=I0TA
      ELSE
       I0TB=I0TA-NTBB*IINTFP
       MXCOR=MXCOR-NTBB*IINTFP
      ENDIF
      CALL GETLST(ICORE(I0TA),1,1,1,1,90)
      IF(IUHF.EQ.1) CALL GETLST(ICORE(I0TB),1,1,1,2,90)
      TAU=.TRUE.
      ENDIF
C
      IF(IUHF.EQ.0) THEN
       IBOT=3
      ELSE
       IBOT=1
      ENDIF
C
      DO 1000 ISPIN=IBOT,3
C
      IF(ISPIN.LT.3) THEN
C
C    AA OR BB CASE.
C
       IF(ISPIN.EQ.1) THEN
        I0T=I0TA
       ELSE
        I0T=I0TB
       ENDIF
C
       LISTT=ISPIN+43
       LISTZ=ISPIN+113
       LISTW=150+ISPIN
       DO 100 IRREP=1,NIRREP
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        I000=1
        I010=I000+DISSYW*NUMSYW*IINTFP
        I020=I010+DISSYT*NUMSYT*IINTFP
        I030=I020+DISSYT*NUMSYT*IINTFP
        IF(MIN(NUMSYT,NUMSYW,DISSYT,DISSYW).NE.0)THEN
         IF(I030.LT.MXCOR)THEN
          CALL V1X2AA(ICORE(I000),ICORE(I010),ICORE(I020),DISSYW,
     &                NUMSYW,DISSYT,NUMSYT,LISTW,LISTT,ITYPE,IRREP,
     &                ICORE(I0T),TAU,POP(1,ISPIN),VRT(1,ISPIN),
     &                LISTZ,ISPIN)
         ELSE
          STOP 'V1X2AA'
         ENDIF
        ENDIF
100    CONTINUE
      ELSE
C
C AB SPIN CASE.  AGAIN, FIRST TRY FOR IN-CORE ALGORITHM.
C
       LISTT=46
       LISTW=153
       DO 110 IRREP=1,NIRREP
C
        DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        I000=1
        I010=I000+NUMSYW*DISSYW*IINTFP
        I020=I010+NUMSYT*DISSYT*IINTFP
        I030=I020+NUMSYT*DISSYT*IINTFP
        IF(MIN(NUMSYT,NUMSYW,DISSYT,DISSYW).NE.0)THEN
         IF(I030.LT.MXCOR)THEN
          CALL V1X2AB(ICORE(I000),ICORE(I010),ICORE(I020),DISSYW,
     &                NUMSYW,DISSYT,NUMSYT,LISTW,LISTT,ITYPE,IRREP,
     &                ICORE(I0TA),ICORE(I0TB),TAU,POP(1,1),POP(1,2),
     &                VRT(1,1),VRT(1,2),IUHF)
C
         ELSE
          CALL INSMEM('V1X2AB',I030,MXCOR)
         ENDIF
        ENDIF
110    CONTINUE
      ENDIF
1000  CONTINUE
      RETURN
      END
