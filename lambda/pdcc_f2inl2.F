      SUBROUTINE PDCC_F2INL2(ICORE,MAXCOR,IUHF,FACT)
C
C THIS ROUTINE CALCULATES THE TERM
C
C  - P(IJ) SUM M L(IM,AB) F(JM)
C
C IN RHF :
C
C - SUM m L(Im,Ab) F(jm) + SUM M L(Mj,Ab) F(JM)
C
C IN UHF
C
C - SUM M L(IM,AB) F(JM) + SUM M L(MJ,AB) F(IM)
C
C - SUM m L(Im,Ab) F(jm) + SUM M L(Mj,Ab) F(IM)
C
C - SUM m L(im,ab) F(jm) + SUM m L(mj,ab) F(im)
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYZ,POP,VRT
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO  
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,
     &             NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /FLAGS/IFLAGS(100)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &              CC2
C
C
      DATA ONE,ONEM,HALF/1.0D0,-1.0D0,0.5D0/
C
C   CALCULATE SIZE OF F(IM) ARRAY
C
      NFAA=NF1AA
      NFBB=NF1BB
      I0AA=MAXCOR+1-NFAA*IINTFP 
      MXCOR=MAXCOR-NFAA*IINTFP
      CALL GETLST(ICORE(I0AA),1,1,1,9,91)
      IF(IUHF.EQ.0) THEN
       I0BB=I0AA
      ELSE
       I0BB=I0AA-NFBB*IINTFP
       MXCOR=MXCOR-NFBB*IINTFP
       CALL GETLST(ICORE(I0BB),1,1,1,10,91)
      ENDIF
C
C     AA AND BB SPIN CASES
C
      IF(IUHF.EQ.1) THEN
C
C
       DO 100 ISPIN=1,2
C
        IF(ISPIN.EQ.1) THEN
         I000=I0AA
        ELSE
         I000=I0BB
        ENDIF
        LISTT=ISPIN+143
        LISTZ=ISPIN+60
C
        DO 50 IRREP=1,NIRREP 
C
        NOCCSQ=0
        DO 45 IRREPJ=1,NIRREP
         IRREPI=DIRPRD(IRREP,IRREPJ)
         NOCCSQ=NOCCSQ+POP(IRREPJ,ISPIN)*POP(IRREPI,ISPIN)
45      CONTINUE
C
         DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT)) 
         DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
         NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
         NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
         I001=1
         I002=I001+IINTFP*NOCCSQ*DISSYT
         I003=I002+IINTFP*NOCCSQ*DISSYZ
         IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
          I004=I003+IINTFP*MAX(DISSYT,DISSYZ)
          IF(I004.LT.MXCOR) THEN
C
C
C    IN CORE VERSION
C
          CALL PDCC_F2L2AA(ICORE(I001),ICORE(I002),ICORE(I000),
     &                POP(1,ISPIN),NOCCSQ,DISSYT,DISSYZ,NUMSYT,
     &                NUMSYZ,NFAA,LISTT,LISTZ,IRREP,ICORE(I003),
     &                FACT)
          ELSE
           STOP 'F2L2AA'
          ENDIF
         ENDIF
50      CONTINUE
100    CONTINUE
      ENDIF
  
C
C      AB SPIN CASE
C
      LISTT=146
      LISTZ=63
C
C   LOOP OVER IRREPS
C
      DO 200 IRREP=1,NIRREP
C
C   RETRIEVE T2 AMPLITUDES AND CALCULATE Z-AMPLITUDES
C
       DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT))
       DISSYZ=IRPDPD(IRREP,ISYTYP(1,LISTZ))
       NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT))
       NUMSYZ=IRPDPD(IRREP,ISYTYP(2,LISTZ))
       I001=1
       I002=I001+IINTFP*NUMSYT*DISSYT
       I003=I002+IINTFP*NUMSYZ*DISSYZ
       IF(MIN(NUMSYT,NUMSYZ,DISSYT,DISSYZ).NE.0) THEN
        I004=I003+IINTFP*MAX(DISSYT,DISSYZ,NUMSYT,NUMSYZ)*3
        IF(I004.LT.MXCOR) THEN
C
C       IN CORE ALGORITHM
C
         CALL PDCC_F2L2AB(ICORE(I001),ICORE(I002),ICORE(I0AA),
     &               ICORE(I0BB),
     &               POP(1,1),POP(1,2),VRT(1,1),VRT(1,2),DISSYT,
     &               DISSYZ,NUMSYT,NUMSYZ,NFAA,NFBB,LISTT,LISTZ,
     &               IRREP,IUHF,ICORE(I003),FACT)
        ELSE
C
C       OUT CORE ALGORITHM
C
         STOP 'F2L2AB'
        ENDIF
       ENDIF
200   CONTINUE
      RETURN
      END
