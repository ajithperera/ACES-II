      SUBROUTINE DFINDT(D,T,F,DISSYD,NUMSYD,DISSYT,
     &                  NUMSYT,NUMR,NUMT,NUMS,IRREPLL,IRREPRR,
     &                  IRREPX,IOFFF,IMOD)
C
C  THIS ROUTINE CALCULATES TERMS OF THE FOLLOWING NATURE :
C
C  A(PQ,RS) = SUM T B(PQ,RT) S(T,S)
C
C  IN-CORE VERSION
C
C  D ...... ARRAY A
C  T ...... ARRAY B
C  F ...... ARRAY S
C  DISSYD . LENGTH OF (P,Q)
C  NUMSYD . LENGTH OF (R,S)
C  DISSYT . LENGTH OF (P,Q)
C  NUMSYT . LENGTH OF (R,T)
C  NUMR ... POPULATION OF R
C  NUMT ... POPULATION OF T
C  NUMS ... POPULATION OF S
C  IRREPLL  IRREP OF A ON THE LEFT SIDE AND OF B ON BOTH SIDES
C  IRREPRR  IRREP OF A ON THE RIGHT HAND SIDE
C  IRREPX . IRREP OF S 
C  IOFFF .. OFFSET VECTOR FOR S
C  IMD .... FLAG IF S (IMOD=1) OR S(TRANSPOSE) (IMOD=2) IS USED
C
CEND
C
C  CODED JAN/91/JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DISSYD,DISSYT,DIRPRD
C
      DIMENSION D(DISSYD,NUMSYD),T(DISSYT,NUMSYT),F(1)
      DIMENSION NUMR(8),NUMT(8),NUMS(8),IOFFF(8),IOFFRR(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA ONE /1.0D0/
C
       IOFFRT=1
C
         IOFFRR(1)=1
         DO 10 IRREP1=1,NIRREP-1
          IRREP2=DIRPRD(IRREP1,IRREPRR)
          IOFFRR(IRREP1+1)=IOFFRR(IRREP1)+NUMS(IRREP1)*NUMR(IRREP2)
10       CONTINUE
C
C  LOOP OVER ALL IRREPS OF T
C
       DO 100 IRREPT=1,NIRREP
C
C  NUMBER OF ORBITALS T
C
        NUMTT=NUMT(IRREPT)
C
C  IRREP OF INDEX R
C
        IRREPR=DIRPRD(IRREPLL,IRREPT)
C
C  NUMBER OF ORBITALS R
C
        NUMRR=NUMR(IRREPR)
C
C  IRREP OF INDEX S
C
        IRREPS=DIRPRD(IRREPT,IRREPX)
C
C  NUMBER OF ORBITALS S
C
        NUMSS=NUMS(IRREPS)
C
C  OFFSET WITHIN R,S
C
        IOFFRS=IOFFRR(IRREPS)
C
C  MULTIPLY F(T,S) WITH T(PQ,RT) ----> (PQR;T) (T,S)
C
        IF(IMOD.EQ.1) THEN
C
C  OFFSET WITHIN F
C
        IFSTART=IOFFF(IRREPS)
C
        CALL XGEMM('N','N',DISSYT*NUMRR,NUMSS,NUMTT,ONE,
     &             T(1,IOFFRT),DISSYT*NUMRR,F(IFSTART),NUMTT,
     &             ONE,D(1,IOFFRS),DISSYT*NUMRR)
        ELSE IF(IMOD.EQ.2) THEN
C
C  OFFSET WITHIN F
C
        IFSTART=IOFFF(IRREPT)
C
        CALL XGEMM('N','T',DISSYT*NUMRR,NUMSS,NUMTT,ONE,
     &             T(1,IOFFRT),DISSYT*NUMRR,F(IFSTART),NUMSS,
     &             ONE,D(1,IOFFRS),DISSYT*NUMRR)
        ELSE
         CALL ERREX
        ENDIF
        
C
C  UPDATE IOFFRS
C
        IOFFRT=IOFFRT+NUMRR*NUMTT
C
100   CONTINUE 
C
      RETURN
      END
