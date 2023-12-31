       SUBROUTINE GAMMA3(ICORE,MAXCOR,IUHF)
C
C THIS ROUTINE COMPUTES THE GAMMA INTERMEDIATE 3:
C
C
C  MBPT(3) : G(IJ,KL) = 1/8 SUM E,F T[1](IJ,EF) T[1](KL,EF)
C
C  MBPT(4) : G(IJ,KL) = 1/8 SUM E,F T2(IJ,EF) T1(KL,EF)
C
C                      + 1/8 SUM E,F T1(IJ,EF) DELTAT2(KL,EF)
C
C WHICH GIVES :
C
C     G(IJ,KL) = 1/16 P(IJ,KL) [2 T2(IJ,EF) - T1(IJ,EF)] T1(KL,EF)
C
C  CCD :     G(IJ,KL) = 1/16 P(IJ,KL) SUM E,F T(IJ,EF) L(KL,EF) 
C
C  QCISD :   G(IJ,KL) = 1/16 P(IJ,KL) SUM E,F T(IJ,EF) L(KL,EF)
C
C  CCSD :    G(IJ,KL) = 1/16 P(IJ,KL) SUM E,F TAU(IJ,EF) L(KL,EF)
C  
C
C  NOTE THAT IN ALL CC METHODS G(IJ,KL) IS CLOSELY RELATED TO
C  V(IJ,KL). ACTUALLY THE RELATIONSHIP
C
C    G(IJ,KL) = 1/8 P(IJ,KL) V(IJ,KL)
C
C  HOLDS. THAT MEANS THAT WE HAVE ONLY TO SYMMETRIZE V(IJ,KL)
C  IN ORDER TO GET G(IJ,KL)
C
C  FOR MBPT(3) AND MBPT(4) HOWEVER, WE HAVE STILL TO CALCULATE
C  G(IJ,KL) SINCE FOR MBPT(3) THERE IS NO V(IJ,KL) INTERMEDIATE
C  AND FPOR MBPT(4) THERE IS A V(IJ,KL) INTERMEDIATE WHICH DIFFERS
C  IN THE AMPLITUDES USED FROM G(IJ,KL)
C
C  NOTE THAT IN THIS ROUINE ACTUALLY 4*G(IJ,KL) IS CALCULATED
C
C  THE FOLLOWING SPIN CASES HAVE TO BE CONSIDERED :
C
C            G(IJ,KL)  T1(IJ,EF) T1(KL,EF)
C
C             AAAA        AAAA      AAAA      (UHF ONLY)
C             BBBB        BBBB      BBBB      (UHF ONLY)
C             ABAB        ABAB      ABAB      (RHF AND UHF)
C
C  THIS ROUTINE USES EXPLICITELY SYMMETRY
C
CEND
C
C CODED AUGUST/90  JG
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,POP,VRT
      LOGICAL ISAME,MBPT4,CC
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
C
      DIMENSION ICORE(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
C
      DATA ONE /1.0D+0/
C
      MXCOR=MAXCOR
C
      CC=CCD.OR.CCSD.OR.QCISD.OR.UCC
      MBPT4=M4DQ.OR.M4SDQ.OR.M4SDTQ
C     
      IF(IUHF.EQ.1) THEN
C
C      AAAA AND BBBB SPIN CASES
C
      DO 1000 ISPIN=1,2   
C
      LISTG=110+ISPIN
      IF(MBPT3) THEN
       LISTT1=43+ISPIN
       LISTT2=43+ISPIN
       ISAME=.TRUE.
       FACT=ONE
      ELSE IF(MBPT4) THEN
       LISTT1=43+ISPIN
       LISTT2=143+ISPIN
       ISAME=.FALSE.
       FACT=ONE
      ELSE IF(CC) THEN
       LISTT1=150+ISPIN
       ISAME=.FALSE.
      ENDIF
C
C LOOP OVER IRREPS.
C
       DO 100 IRREP=1,NIRREP
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT1))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT1))
        I001=1
        IF(.NOT.CC) THEN
        IF(ISAME) THEN
         I002=I001
        ELSE
         I002=I001+IINTFP*NUMSYT*DISSYT
        ENDIF
        I003=I002+IINTFP*NUMSYT*DISSYT
        ELSE
         I002=I001
         I003=I002
        ENDIF
        IF(MIN(NUMSYT,DISSYT).NE.0)THEN
         I004=I003+IINTFP*NUMSYT*NUMSYT
         IF(I004.LT.MXCOR) THEN
C  
C         IN CORE VERSION
C
          CALL G3ALL(ICORE(I001),ICORE(I002),ICORE(I003),
     &               CC,MBPT4,FACT,ISPIN,ISAME,DISSYT,
     &               NUMSYT,LISTT1,LISTT2,LISTG,IRREP)
         ELSE
C
C         OUT OF CORE VERSION
C
          CALL INSMEM('G3ALL',I004,MXCOR)
         ENDIF
        ENDIF
100    CONTINUE
1000   CONTINUE
      ENDIF
C
C       AB SPIN CASE
C
      LISTG=113
      IF(MBPT3) THEN
       LISTT1=46  
       LISTT2=46
       ISAME=.TRUE.
       FACT=ONE
      ELSE IF(MBPT4) THEN
       LISTT1=46
       LISTT2=146
       ISAME=.FALSE.
       FACT=ONE
      ELSE IF(CC) THEN
       LISTT1=153
       ISAME=.FALSE.
      ENDIF
C
C      LOOP OVER IRREPS.
C
       DO 110 IRREP=1,NIRREP
C
        DISSYT=IRPDPD(IRREP,ISYTYP(1,LISTT1))
        NUMSYT=IRPDPD(IRREP,ISYTYP(2,LISTT1))
        I001=1
        IF(.NOT.CC) THEN
         IF(ISAME) THEN
          I002=I001
         ELSE
          I002=I001+IINTFP*NUMSYT*DISSYT
         ENDIF
         I003=I002+IINTFP*NUMSYT*DISSYT
        ELSE
         I002=I001
         I003=I002
        ENDIF
        IF(MIN(NUMSYT,DISSYT).NE.0)THEN
         I004=I003+IINTFP*NUMSYT*NUMSYT
         IF(I004.LT.MXCOR) THEN
C
C         IN CORE VERSION
C
          CALL G3ALL(ICORE(I001),ICORE(I002),ICORE(I003),CC,MBPT4,
     &               FACT,3,ISAME,DISSYT,NUMSYT,LISTT1,LISTT2,
     &               LISTG,IRREP)
         ELSE
          CALL INSMEM('G3ALL',I004,MXCOR) 
         ENDIF
        ENDIF
110    CONTINUE
C
      if(iuhf.eq.0) then
      call checkgam1(icore,13,113,ONE,IUHF,1,POP)
      endif
      if(iuhf.eq.1) then
       call checkgam(icore,13,113,ONE)
       call checkgam(icore,11,111,ONE)
       call checkgam(icore,12,112,ONE)
      endif

      RETURN
      END
