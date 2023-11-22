      SUBROUTINE DFINDT3(D,T,F,ISTART,IEND,DISSYD,NUMSYD,DISSYT,
     &                   NUMSYT,NUMR,NUMT,NUMS,IRREPLL,IRREPRR,
     &                   IRREPX,IOFFF,IMOD)
C
C  THIS ROUTINE CALCULATES TERMS OF THE FOLLOWING NATURE :
C
C  A(PQ,SR) = SUM T B(PQ,TR) S(T,S)
C
C  OUT-OF-CORE VERSION USING PARTIAL T-ARRAY
C  THIS VERSION DEALS WITH TRUNCATED RIGHT-HAND SIDE.
C
C  D ...... ARRAY A
C  T ...... ARRAY B
C  F ...... ARRAY S
C  ISTART.. FIRST DISTRIBUTION OF D
C  IEND.... LAST  DISTRIBUTION OF D
C  DISSYD . FULL LENGTH OF (P,Q)
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
      DIMENSION NUMR(8),NUMT(8),NUMS(8),IOFFF(8)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
C
      DATA ONE /1.0D0/
C
      IOFFRT=1
      IOFFRS=1
C
C  LOOP OVER ALL IRREPS OF T
C
      DO 10 IRREPR=1,NIRREP
C
C  NUMBER OF ORBITALS R
C
       NUMRR=NUMR(IRREPR)
C
C  IRREP OF INDEX T
C
       IRREPT=DIRPRD(IRREPLL,IRREPR)
C
C  NUMBER OF ORBITALS T
C
       NUMTT=NUMT(IRREPT)
C
C  IRREP OF INDEX S
C
       IRREPS=DIRPRD(IRREPT,IRREPX)
C
C  NUMBER OF ORBITALS S
C
       NUMSS=NUMS(IRREPS)
C
       IRSTART=IOFFRT
       IREND=IOFFRT+NUMTT*NUMRR-1
       IF(ISTART.LE.IREND.AND.IEND.GE.IRSTART) THEN
C
C  MULTIPLY F(T,S) WITH T(PQ,RT) ----> (PQR;T) (T,S)
C
        IF(IMOD.EQ.1) THEN
C
C  OFFSET WITHIN F
C
         IFSTART=IOFFF(IRREPS)
C
         DO 100 IR=1,NUMRR
C
          IRSTART=IOFFRT+(IR-1)*NUMTT
          IREND=IOFFRT+IR*NUMTT-1
          IF(ISTART.LE.IREND.AND.IEND.GE.IRSTART) THEN
C
           ISTART1=MAX(ISTART,IRSTART)
           IEND1=MIN(IEND,IREND) 
           NUMT1=IEND1-ISTART1+1
           IT=ISTART1-IRSTART
           IFSTART1=IFSTART+IT
C
           CALL XGEMM('N','N',DISSYT,NUMSS,NUMT1,ONE,
     &                 T(1,IOFFRT+(IR-1)*NUMTT+IT-ISTART+1),
     &                 DISSYT,F(IFSTART1),NUMTT,
     &                 ONE,D(1,IOFFRS+(IR-1)*NUMSS),
     &                 DISSYD)
C
          ENDIF
C
100      CONTINUE
C
        ELSE IF(IMOD.EQ.2) THEN
C
C  OFFSET WITHIN F
C
         IFSTART=IOFFF(IRREPT)
C
         DO 200 IR=1,NUMRR 
C
          IRSTART=IOFFRT+(IR-1)*NUMTT
          IREND=IOFFRT+IR*NUMTT-1
          IF(ISTART.LE.IREND.AND.IEND.GE.IRSTART) THEN
C
           ISTART1=MAX(ISTART,IRSTART)
           IEND1=MIN(IEND,IREND) 
           NUMT1=IEND1-ISTART1+1
           IT=ISTART1-IRSTART
           IFSTART1=IFSTART+IT*NUMSS
C
           CALL XGEMM('N','T',DISSYT,NUMSS,NUMT1,ONE,
     &                T(1,IOFFRT+(IR-1)*NUMTT+IT-ISTART+1),
     &                DISSYT,F(IFSTART1),NUMSS,
     &                ONE,D(1,IOFFRS+(IR-1)*NUMSS),
     &                DISSYD)
C
          ENDIF
C
200      CONTINUE
C
        ELSE
C
         CALL ERREX
C
        ENDIF
C
       ENDIF        
C
C  UPDATE IOFFRS AND IOFFRT
C
       IOFFRT=IOFFRT+NUMRR*NUMTT
       IOFFRS=IOFFRS+NUMRR*NUMSS
C
10    CONTINUE 
C
      RETURN
      END