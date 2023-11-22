      SUBROUTINE FORMU(IRREP,NPERT,N,UAI,EVAL,POP,VRT,NOCC)
C
C   THIS ROUTINE CALCULATES 
C
C   U(AI) = U(AI)/(EI-EA)
C
C REQUIRED FOR SOLVING THE CPHF EQUATION
C
C THIS ROUTINE USES EXPLICITELY SYMMETRY
C 
CEND
C
C CODED AUGUST/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      DIMENSION UAI(N,NPERT),EVAL(1),POP(8),VRT(8),IV(8)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
C  CALCULATE FIRST OFFSET FOR THE VIRTUAL ORBITALS
C
      IV(1)=NOCC
      DO 1 IRREPJ=1,NIRREP-1
      IV(IRREPJ+1)=IV(IRREPJ)+VRT(IRREPJ)
1     CONTINUE
C
      IND0=0
      IOFFP=0
C
C  LOOP OVER THE IRREP OF THE OCCUPIED ORBITALS, THE
C  IRREP OF THE VIRTUAL ORBITALS IS THEN SIMPLY GIVEN
C  AS THE DIRECT PRODUCT OF IRREP WITh IRREPJ
C
      DO 10 IRREPJ=1,NIRREP
       IRREPI=DIRPRD(IRREP,IRREPJ)
       NOCCI=POP(IRREPJ)
       NVRTI=VRT(IRREPI)
       IOFFV=IV(IRREPI)
       DO 6 I=1,NOCCI
        INDI=I+IOFFP
        DO 5 IPERT=1,NPERT
        IND=IND0
        DO 5 IA=1,NVRTI
        INDA=IA+IOFFV
        IND=IND+1
        UAI(IND,IPERT)=UAI(IND,IPERT)/(EVAL(INDI)-EVAL(INDA))
5      CONTINUE 
       IND0=IND0+NVRTI
6      CONTINUE
       IOFFP=IOFFP+NOCCI
10    CONTINUE
      RETURN
      END