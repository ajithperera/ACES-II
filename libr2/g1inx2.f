      SUBROUTINE G1INX2(ICORE,MAXCOR,IUHF,IOFF)
C
C THIS ROUTINE CALCULATES THE TERM
C
C  P(AB) SUM E TAU(IJ,AE) G(E,B)
C
C USING SYMMETRY PACKED ARRAYS
C 
C
C IN RHF :
C
C SUM e TAU(Ij,Ae) G(be) - SUM E  TAU(Ij,Eb) G(AE)
C
C IN UHF
C
C SUM E TAU(IJ,AE) G(BE) - SUM E TAU(IJ,EB) G(AE)
C
C SUM e TAU(Ij,Ae) G(be) - SUM E TAU(Ij,Eb) G(AE)
C
C SUM e TAU(ij,ae) G(be) - SUM e TAU(ij,eb) G(ae)
C
C 
C FOR QCISD AND CCSD IN ADDITION THE FOLLWING TERM
C MUST BE ADDED TO G(AE)
C
C  G(AE) = G(AE) - SUM M L(M,E) T(M,A)
C
CEND
C
C
C CODED JG AUGUST/90
C
C
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
      DATA ONE,ONEM /1.0D0,-1.0D0/
C
      TAU=CCSD
C
C    CALCULATE SIZE OF G(A,E) ARRAY
C    AND GET THESE ARRAYS
C
      NFAA=NF2AA
      NFBB=NF2BB
      I0AA=MAXCOR+1-NFAA*IINTFP
      MXCOR=MAXCOR-NFAA*IINTFP
      CALL GETLST(ICORE(I0AA),1,1,1,1+IOFF,192)
      IF(IUHF.EQ.0) THEN
       I0BB=I0AA
      ELSE
       I0BB=I0AA-NFBB*IINTFP
       MXCOR=MXCOR-NFBB*IINTFP
       CALL GETLST(ICORE(I0BB),1,1,1,2+IOFF,192)
      ENDIF
C
C FOR QCISD AND CCSD METHODS GET THE T1 AND L1 AMPLITUDES
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
C FORM ADDITIONAL CONTRIBUTION TO H(EA)
C
       DO 300 ISPIN=1,IUHF+1
C
C  GET OFFSETS FOR T1, L1, AND H1
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
        CALL XGEMM('N','T',NVRT,NVRT,NOCC,ONEM,ICORE(IOFFT),NVRT,
     &             ICORE(IOFFL),NVRT,ONE,ICORE(IOFFF),NVRT)
        IOFFF=IOFFF+IINTFP*NVRT*NVRT
        IOFFT=IOFFT+IINTFP*NOCC*NVRT
        IOFFL=IOFFL+IINTFP*NOCC*NVRT
20     CONTINUE
300    CONTINUE
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
         CALL G1X2AA(ICORE(I001),ICORE(I002),ICORE(I002),
     &               ICORE(I001),ICORE(I0T),TAU,ISPIN,ICORE(I000),
     &               POP(1,ISPIN),VRT(1,ISPIN),NVRTSQ,DISSYT,
     &               DISSYZ,NUMSYT,NUMSYZ,NFAA,LISTT,LISTZ,IRREP,
     &               ICORE(I003))
        ELSE
         STOP 'G1X2AA'
        ENDIF
       ENDIF
50    CONTINUE
100   CONTINUE
      ENDIF
  
C
C      AB SPIN CASE
C
       LISTT=46
       LISTZ=116
C
C    LOOP OVER IRREPS
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
         I004=I003+IINTFP*MAX(NUMSYT,NUMSYZ,DISSYZ,DISSYT)*3
         IF(I004.LT.MXCOR) THEN
C
C       IN CORE ALGORITHM
C
         CALL G1X2AB(ICORE(I001),ICORE(I002),ICORE(I002),ICORE(I001),
     &               ICORE(I0TA),ICORE(I0TB),TAU,ICORE(I0AA),
     &               ICORE(I0BB),POP(1,1),POP(1,2),VRT(1,1),
     &               VRT(1,2),DISSYT,DISSYZ,NUMSYT,NUMSYZ,NFAA,NFBB,
     &               LISTT,LISTZ,IRREP,IUHF,ICORE(I003))
        ELSE
C
C       OUT CORE ALGORITHM
C
        STOP 'G1X2AB'
        ENDIF
       ENDIF
200   CONTINUE
      RETURN
      END
