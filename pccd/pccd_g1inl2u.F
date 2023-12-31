      SUBROUTINE PCCD_G1INL2U(ICORE,MAXCOR,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYZ,POP,VRT
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      DATA ONE,HALFM /1.0D0,-0.5D0/
C
      NFAA=NF2AA
      NFBB=NF2BB
      I0AA=MAXCOR+1-NFAA*IINTFP
      MXCOR=MAXCOR-NFAA*IINTFP

      CALL GETLST(ICORE(I0AA),1,1,1,1,192)

      IF(IUHF.EQ.0) THEN
       I0BB=I0AA
      ELSE
       I0BB=I0AA-NFBB*IINTFP
       MXCOR=MXCOR-NFBB*IINTFP
       CALL GETLST(ICORE(I0BB),1,1,1,2,192) 
C
C AUGMENT G IF THIS IS CCSD OR QCISD
C
      ENDIF
C
C     AA AND BB SPIN CASES
C
      IF(IUHF.EQ.1) THEN
C
C      THESE CASES ARE ONLY NECCESARY IN THE UHF CASE
C      IN RHF THE AAAA AMPLITUDES ARE CALCULATED FROM
C      THE ABAB AMPLITUDES
C
       DO 100 ISPIN=1,2
C
        IF(ISPIN.EQ.1) THEN 
         I000=I0AA
         NFSIZE=NFAA
        ELSE
         I000=I0BB
         NFSIZE=NFBB
        ENDIF
        LISTT=ISPIN+13
        LISTZ=ISPIN+60
C
        DO 50 IRREP=1,NIRREP
C
        NVRTSQ=0
        DO 45 IRREPJ=1,NIRREP
         IRREPI=DIRPRD(IRREP,IRREPJ)
         NVRTSQ=NVRTSQ+VRT(IRREPJ,ISPIN)*VRT(IRREPI,ISPIN)
45      CONTINUE
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
        I001=1
        I002=I001+IINTFP*NUMSYT*NVRTSQ
        I003=I002+IINTFP*NUMSYZ*NVRTSQ
        IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
         I004=I003+IINTFP*MAX(NUMSYT,NUMSYZ)
         IF(I004.LT.MXCOR) THEN
C
C     IN CORE VERSION
C
         CALL G1L2AA(ICORE(I001),ICORE(I002),ICORE(I002),
     &               ICORE(I001),ICORE(I000),VRT(1,ISPIN),
     &               NVRTSQ,DISSYT,DISSYZ,NUMSYT,
     &               NUMSYZ,NFSIZE,LISTT,LISTZ,IRREP,ICORE(I003))
        ELSE
        CALL INSMEM("G1INL2U",I004,MAXCOR)
        ENDIF
       ENDIF
50    CONTINUE
100   CONTINUE
      ENDIF
  
C
C      AB SPIN CASE
C
       LISTT=16
       LISTZ=63
C
C    LOOP OVER IRREPS
C
       DO 200 IRREP=1,NIRREP
C
C    RETRIEVE AMPLITUDES AND CALCULATE CONTRIBUTION TO Z
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
        DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
        NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
        I001=1
        I002=I001+IINTFP*NUMSYT*DISSYT
        I003=I002+IINTFP*NUMSYZ*DISSYZ
        IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
         I004=I003+IINTFP*MAX(NUMSYT,NUMSYZ,DISSYZ,DISSYT)*3
         IF(I004.LT.MXCOR) THEN
C
C       IN CORE ALGORITHM
C
         CALL G1L2AB(ICORE(I001),ICORE(I002),ICORE(I002),ICORE(I001),
     &               ICORE(I0AA),ICORE(I0BB),
     &               POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),DISSYT,
     &               DISSYZ,NUMSYT,NUMSYZ,NFAA,NFBB,
     &               LISTT,LISTZ,IRREP,IUHF,ICORE(I003))
        ELSE
C
C       OUT CORE ALGORITHM
C
        CALL INSMEM("G1INL2U",I004,MAXCOR)
        ENDIF
       ENDIF
200   CONTINUE
      RETURN
      END
