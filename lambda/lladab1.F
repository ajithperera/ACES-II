      SUBROUTINE LLADAB1(W,T,Z,MAXSIZE,DISSYW,NUMSYW,DISSYT,
     &                   NUMSYT,LISTW,LISTT,ITYPE,IRREP,CC,CCSD)
C
C THIS ROUTINE FORMS THE AB LADDER TYPE CONTRACTION FOR MBPT(3) OR 
C THE GENERALIZED LADDER CONTRACTIONS IN THE INTERMEDIATE-DRIVEN
C COUPLED-CLUSTER SCHEME.  THIS CURRENTLY RUNS ONLY FOR UHF CALCULATIONS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYW,DISREAD,DISLEFT,DISMAX,DIRPRD
      INTEGER POP,VRT
      LOGICAL CC,CCSD
C
      DIMENSION W(DISSYW,NUMSYW),T(DISSYT,NUMSYT)
      DIMENSION Z(DISSYT,NUMSYT)
C
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/FLAGS/IFLAGS(100)
C
      DATA ZILCH,ONE /0.D0,1.0D0/
C
C PICK UP THE RELEVANT T2 AND W PIECES.
C
      CALL GETLST(T,1,NUMSYT,1,IRREP,LISTT)

C
C CAN WE DO IT IN CORE ?
C
      IF(MAXSIZE.GE.DISSYW*NUMSYW) THEN
C
C IN CORE
C
      CALL GETLST(W,1,NUMSYW,2,IRREP,LISTW)
C
C COMPUTE THE LADDER CONTRIBUTION.
C
      IF(ITYPE.EQ.6)THEN
C
       CALL XGEMM('T','N',DISSYW,NUMSYT,DISSYW,ONE,W,
     &             DISSYW,T,DISSYT,ZILCH,Z,DISSYT)
      ELSE
C
C MULTIPLY WITH THE TRANSPOSE OF THE INTERMEDIATE W(MN,IJ)
C
       CALL XGEMM('N','T',DISSYT,DISSYW,DISSYW,ONE,
     &            T,DISSYT,W,DISSYW,ZILCH,Z,DISSYT)
C
C FOR CC METHODS MULTIPLY ALSO <MN||AB> WITH V(MN,IJ)
C
       IF(CC)THEN
        CALL GETLST(T,1,NUMSYT,1,IRREP,16)
        CALL GETLST(W,1,NUMSYW,2,IRREP,153)
        CALL XGEMM('N','N',DISSYT,NUMSYW,DISSYW,ONE,
     &             T,DISSYT,W,DISSYW,ONE,Z,DISSYT)
       ENDIF
      ENDIF
C
      ELSE
C
C  WE HAVE TO DO IT OUT OF CORE
C
      CALL ZERO(Z,NUMSYT*DISSYT)
C
      DISMAX=MAXSIZE/DISSYW
      IOFFSET=1
      DISLEFT=NUMSYW
10    CONTINUE
       DISREAD=MIN(DISMAX,DISLEFT)
       DISLEFT=DISLEFT-DISREAD
       CALL GETLST(W,IOFFSET,DISREAD,2,IRREP,LISTW)
C
C  COMPUTE NOW THE LADDER CONTRIBUTION
C 
       IF(ITYPE.EQ.6) THEN
        CALL XGEMM('T','N',DISREAD,NUMSYT,DISSYW,ONE,W,DISSYW,
     &             T,DISSYT,ONE,Z(IOFFSET,1),DISSYT)
       ELSE
        CALL XGEMM('N','T',DISSYT,DISSYW,DISREAD,ONE,T(1,IOFFSET),
     &             DISSYT,W,DISSYW,ONE,Z,DISSYT)
       ENDIF
       IOFFSET=IOFFSET+DISREAD
      IF(DISLEFT.NE.0) GO TO 10
C
C    IF ITYPE.EQ.1 CALCULATE ADDITIONAL CONTRIBUTION FOR CC
C
      IF(ITYPE.EQ.1.AND.CC) THEN
       CALL GETLST(T,1,NUMSYT,1,IRREP,16)
       IOFFSET=1
       DISLEFT=NUMSYW
20     CONTINUE
        DISREAD=MIN(DISMAX,DISLEFT)
        DISLEFT=DISLEFT-DISREAD
        CALL GETLST(W,IOFFSET,DISREAD,2,IRREP,153)
        CALL XGEMM('N','N',DISSYT,DISREAD,DISSYW,ONE,
     &             T,DISSYT,W,DISSYW,ONE,Z(1,IOFFSET),DISSYT)
        IOFFSET=IOFFSET+DISREAD
       IF(DISLEFT.NE.0) GO TO 20
      ENDIF
      ENDIF
C
C CALCULATE CONTRIBUTION TO L1 INCREMENT FOR CC CASES (ARISES BECAUSE
C THE ABCD CONTRIBUTION TO HBAR(AIBC) IS NOT CALCULATED ANY MORE
C
c     IF(CC.AND.ITYPE.EQ.6.and.idoit.eq.33)THEN
      IF(CCSD.AND.ITYPE.EQ.6)THEN
       IF(IFLAGS(3).NE.2) THEN
        CALL SCOPY(DISSYT*NUMSYT,Z,1,T,1)
        I0T=1
        ITAR=I0T+MAX(NT(1),NT(2))
        ISCR=ITAR+MAX(NT(1),NT(2))
        CALL GETLST(W(I0T,1),1,1,1,2,90)
        CALL GETLST(W(ITAR,1),1,1,1,3,90)
        CALL DOT24(IRREP,W(ITAR,1),W(I0T,1),T,
     &             W(ISCR,1),DISSYT,VRT(1,1),POP(1,1),
     &             VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),
     &             'TSTS')
        CALL PUTLST(W(ITAR,1),1,1,1,3,90)
        CALL GETLST(W(I0T,1),1,1,1,1,90)
        CALL GETLST(W(ITAR,1),1,1,1,4,90)
        CALL DOT24(IRREP,W(ITAR,1),W(I0T,1),T,
     &             W(ISCR,1),DISSYT,VRT(1,2),POP(1,2),
     &             VRT(1,1),VRT(1,2),POP(1,1),POP(1,2),
     &             'STST')
        CALL PUTLST(W(ITAR,1),1,1,1,4,90)
       ENDIF
      ENDIF
C
C AUGMENT T2 INCREMENT.
C
      CALL SUMSYM2(Z,T,NUMSYT*DISSYT,1,IRREP,63)
C
      RETURN
      END
