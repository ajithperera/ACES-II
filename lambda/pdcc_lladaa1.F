      SUBROUTINE PDCC_LLADAA1(W,T,Z,MAXSIZE,DISSYW,NUMSYW,DISSYT,
     &                   NUMSYT,LISTW,LISTT,ITYPE,IRREP,LISTZ,
     &                   ISPIN,CC,CCSD)
C
C THIS ROUTINE COMPUTES THE A-A OR B-B LADDER CONTRIBUTIONS TO THE
C L2 INCREMENT WHEN SOLVING THE LAMBDA EQUATIONS
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYW,DISLEFT,DISREAD,DISMAX,DIRPRD
      INTEGER DISSYF,NUMDSF,POP,VRT
      LOGICAL CC,CCSD
C
      DIMENSION W(DISSYW,NUMSYW),T(DISSYT,NUMSYT)
      DIMENSION Z(DISSYT,NUMSYT)
C
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/FLAGS/IFLAGS(100)
C
      DATA ZILCH,ONE /0.0D0,1.0D0/
C
C PICK UP THE RELEVANT T2 AND W PIECES.
C
      CALL GETLST(T,1,NUMSYT,1,IRREP,LISTT)
C
C  CAN WE DO IT INCORE ?
C
      IF(MAXSIZE.GE.DISSYW*NUMSYW) THEN
C
C  INCORE
C
      CALL GETLST(W,1,NUMSYW,2,IRREP,LISTW)
C
C COMPUTE THE LADDER CONTRIBUTION.
C
      IF(ITYPE.EQ.6)THEN
        CALL XGEMM('T','N',DISSYW,NUMSYT,DISSYT,ONE,
     &             W,DISSYW,T,DISSYT,ZILCH,Z,DISSYW)
      ELSE
C
C  MULTIPLY L2 AMPLITUDES WITH THE TRANSPOSE OF W(MN,IJ)
C
       CALL XGEMM('N','T',DISSYT,DISSYW,DISSYW,ONE,
     &            T,DISSYT,W,DISSYW,ZILCH,Z,DISSYT)
C
C  FOR CC METHODS MULTIPLY IN ADDITION <MN//AB> WITH V(MN,IJ)
C
       IF(CC)THEN
        CALL GETLST(T,1,NUMSYT,1,IRREP,13+ISPIN)
        CALL GETLST(W,1,NUMSYW,2,IRREP,150+ISPIN)
        CALL XGEMM('N','N',DISSYT,NUMSYW,DISSYW,ONE,
     &             T,DISSYT,W,DISSYW,ONE,Z,DISSYT)
      ENDIF
      ENDIF
C
      ELSE
C
C   WE HAVE TO DO IT OUT OF CORE
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
C   COMPUTE NOW THE LADDER CONTRIBUTION
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
C CALCULATE ADDITIONAL CONTRIBUTION FOR HHL CALCULATION    
C
      IF(ITYPE.EQ.1.AND.CC) THEN
       CALL GETLST(T,1,NUMSYT,1,IRREP,13+ISPIN)     
       IOFFSET=1
       DISLEFT=NUMSYW
20     CONTINUE
        DISREAD=MIN(DISMAX,DISLEFT)
        DISLEFT=DISLEFT-DISREAD
        CALL GETLST(W,IOFFSET,DISREAD,2,IRREP,150+ISPIN)
        CALL XGEMM('N','N',DISSYT,DISREAD,DISSYW,ONE,
     &             T,DISSYT,W,DISSYW,ONE,Z(1,IOFFSET),DISSYT)
        IOFFSET=IOFFSET+DISREAD
        IF(DISLEFT.NE.0) GO TO 20
      ENDIF
      ENDIF
C
C FORM CONTRIBUTION TO L1 INCREMENT FOR CC PPL CASE
C
      IF(CCSD.AND.ITYPE.EQ.6)THEN
       IF(IFLAGS(3).NE.2) THEN
        DISSYF=IRPDPD(IRREP,18+ISPIN)
        NUMDSF=IRPDPD(IRREP,20+ISPIN)
        CALL SCOPY(NUMSYT*DISSYT,Z,1,W,1)
        CALL SYMEXP(IRREP,POP(1,ISPIN),DISSYT,W)
        CALL SYMEXP2(IRREP,VRT(1,ISPIN),DISSYF,DISSYT,NUMDSF,T,W)
        I0T=1
        ITAR=I0T+NT(ISPIN)
        ISCR=ITAR+NT(ISPIN)
        CALL GETLST(W(I0T,1),1,1,1,ISPIN,90)   
        CALL GETLST(W(ITAR,1),1,1,1,2+ISPIN,90)   
        CALL DOT24(IRREP,W(ITAR,1),W(I0T,1),T,W(ISCR,1),
     &             DISSYF,VRT(1,ISPIN),POP(1,ISPIN),
     &             VRT(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),POP(1,ISPIN),
     &             'TSTS')
        CALL PUTLST(W(ITAR,1),1,1,1,2+ISPIN,90)
       ENDIF
      ENDIF
C
C AUGMENT T2 INCREMENT.
C
      CALL SUMSYM2(Z,T,NUMSYT*DISSYT,1,IRREP,LISTZ)
C
      RETURN
      END
