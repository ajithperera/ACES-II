      SUBROUTINE G2INX2(ICORE,MAXCOR,IUHF,IOFF)
C
C THIS PROGRAM CALCULATES THE TERM
C
C -  P(IJ) SUM M TAU(IM,AB) G(JM)
C
C SYMMETRY ADAPTED 
C 
C
C IN RHF :
C
C - SUM m TAU(Im,Ab) G(jm) + SUM M TAU(Mj,Ab) G(IM)
C
C IN UHF
C
C - SUM M TAU(IM,AB) G(JM) + SUM M TAU(MJ,AB) G(IM)
C
C - SUM m TAU(Im,Ab) G(jm) + SUM M TAU(Mj,Ab) G(IM)
C
C - SUM m TAU(im,ab) G(jm) + SUM m TAU(mj,ab) G(im)
C
C
C  FOR CCSD AND QCISD THE FOLLOWING TERM HAS TO BE ADDED TO G(IM)
C 
C  G(JM) = G(JM) + SUM E L(M,E) T(J,E)
C
CEND
C
C CODED JG JUNE/90
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,DISSYZ,POP,VRT
      LOGICAL TAU
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      DIMENSION ICORE(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWO  
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
C
      DATA ONE,ONEM/1.0D0,-1.0D0/
C
      TAU=CCSD
C
C   CALCULATE SIZE OF F(M,I) ARRAY
C
      NFAA=NF1AA
      NFBB=NF1BB
      I0AA=MAXCOR+1-NFAA*IINTFP 
      MXCOR=MAXCOR-NFAA*IINTFP
      CALL GETLST(ICORE(I0AA),1,1,1,1+IOFF,191)
      IF(IUHF.EQ.0) THEN
       I0BB=I0AA
      ELSE
       I0BB=I0AA-NFBB*IINTFP
       MXCOR=MXCOR-NFBB*IINTFP
       CALL GETLST(ICORE(I0BB),1,1,1,2+IOFF,191)
      ENDIF
C
C  FOR QCISD AND CCSD METHODS GET THE T1 AND L1 AMPLITUDES
C
      IF(QCISD.OR.CCSD) THEN
       I0TA=I0BB-NTAA*IINTFP
       I0LA=I0TA-NTAA*IINTFP
       MXCOR=MXCOR-2*NTAA*IINTFP
       CALL GETLST(ICORE(I0TA),1,1,1,1,90)
       CALL GETLST(ICORE(I0LA),1,1,2,1,190)
       IF(IUHF.EQ.0) THEN
        I0TB=I0TA
        I0LB=I0LA
       ELSE
        I0TB=I0LA-NTBB*IINTFP
        I0LB=I0TB-NTBB*IINTFP
        MXCOR=MXCOR-2*NTBB*IINTFP
        CALL GETLST(ICORE(I0TB),1,1,1,2,90)
        CALL GETLST(ICORE(I0LB),1,1,2,2,190)
       ENDIF
C
C FORM ADDITIONAL CONTRIBUTION TO H(MI)
C
       DO 300 ISPIN=1,IUHF+1
C
C  GET OFFSETS FOR T1, L1 AND H1
C
        IF(ISPIN.EQ.1) THEN
         IOFFT=I0TA
         IOFFL=I0LA
         IOFFF=I0AA
        ELSE
         IOFFT=I0TB
         IOFFL=I0LB
         IOFFF=I0BB
        ENDIF
C
C  LOOP OVER IRREPS
C
        DO 20 IRREP=1,NIRREP
         NOCC=POP(IRREP,ISPIN)
         NVRT=VRT(IRREP,ISPIN)
         CALL XGEMM('T','N',NOCC,NOCC,NVRT,ONE,ICORE(IOFFT),NVRT,
     &              ICORE(IOFFL),NVRT,ONE,ICORE(IOFFF),NOCC)
         IOFFF=IOFFF+IINTFP*NOCC*NOCC
         IOFFT=IOFFT+IINTFP*NOCC*NVRT
         IOFFL=IOFFL+IINTFP*NOCC*NVRT
20      CONTINUE
300     CONTINUE
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
         I0T=I0TA
        ELSE
         I000=I0BB
         I0T=I0TB
        ENDIF
        LISTT=ISPIN+43
        LISTZ=ISPIN+113
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
C    IN CORE VERSION
C
          CALL G2X2AA(ICORE(I001),ICORE(I002),ICORE(I0T),TAU,ISPIN,
     &                ICORE(I000),POP(1,ISPIN),VRT(1,ISPIN),NOCCSQ,
     &                DISSYT,DISSYZ,NUMSYT,NUMSYZ,NFAA,LISTT,LISTZ,
     &                IRREP,ICORE(I003))
          ELSE
           STOP 'G2X2AA'
          ENDIF
         ENDIF
50      CONTINUE
100    CONTINUE
      ENDIF
  
C
C      AB SPIN CASE
C
      LISTT=46
      LISTZ=116
C
C   LOOP OVER IRREPS
C
      DO 200 IRREP=1,NIRREP
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
         CALL G2X2AB(ICORE(I001),ICORE(I002),ICORE(I0TA),ICORE(I0TB),
     &               TAU,ICORE(I0AA),ICORE(I0BB),POP(1,1),POP(1,2),
     &               VRT(1,1),VRT(1,2),DISSYT,DISSYZ,NUMSYT,NUMSYZ,
     &               NFAA,NFBB,LISTT,LISTZ,IRREP,IUHF,ICORE(I003))
        ELSE
C
C       OUT CORE ALGORITHM
C
         STOP 'G2X2AB'
        ENDIF
       ENDIF
200   CONTINUE
      RETURN
      END
