      SUBROUTINE GETD(DOOA,DOOB,DVVA,DVVB,DVOA,DVOB,IUHF,ROHF)
C
C THIS SUBROUTINE SAVES THE RELAXED DENSITY MATRIX ON THE
C GAMLAM FILE. THIS IS ONLY REQUIRED FOR SECOND DERIVATIVE
C CALCULATIONS.
C
C THE LISTS READ ARE:
C
C  DOOA    1,160
C  DOOB    2,160    UHF AND ROHF ONLY
C  DVVA    3,160
C  DVVB    4,160    UHF AND ROHF ONLY
C  DVOA    5,160  
C  DVOB    6,160    UHF AND ROHF ONLY
C
CEND
C
C CODED APRIL/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL ROHF 
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      INTEGER POP,VRT
      DIMENSION DOOA(*),DOOB(*),DVVA(*),DVVB(*),DVOA(*),DVOB(*)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
C
      CALL GETLST(DOOA,1,1,1,1,160)
      CALL GETLST(DVVA,1,1,1,3,160)
      IF(ROHF) CALL GETLST(DVOA,1,1,1,5,160)
C
      IF(IUHF.NE.0) THEN
C
       CALL GETLST(DOOB,1,1,1,2,160)
       CALL GETLST(DVVB,1,1,1,4,160)
       IF(ROHF) CALL GETLST(DVOB,1,1,1,6,160)
C
      ENDIF
C
      RETURN
      END