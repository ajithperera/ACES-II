      SUBROUTINE FNDLEFT(LENGTH,TEST,VINPUT,VEXCLUDE,TMP,NEXCLUDE,
     &                   VLOC,ILOC,TOL,THRESH)
C
C FINDS MINIMUM ELEMENT IN VECTOR VINPUT WHICH IS NOT WITHIN TOL
C OF THE VALUES VEXCLUDE(1..NEXCLUDE), AND WHICH ARE LARGER THAN
C THRESH (IN CASE OF CORE-EXCITATIONS).
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VINPUT(LENGTH),VEXCLUDE(NEXCLUDE),TMP(LENGTH)
C
      CALL SCOPY(LENGTH,VINPUT,1,TMP,1)
      DO 5 I = 1, LENGTH
         IF (TMP(I).LT. THRESH) TMP(I) = 1.D+30
 5    CONTINUE
      DO 10 I=1,NEXCLUDE
       CALL FNDCLOSE(LENGTH,TMP,VEXCLUDE(I),X,ILOC)
       X2=ABS(VEXCLUDE(I)-TMP(ILOC))
       IF(X2.LT.TOL)TMP(ILOC)=1.D+30
10    CONTINUE
C
      XCLOSE=1.D+30
      DO 20 I=1,LENGTH
       X=ABS(TMP(I)-TEST)
       IF(X.LT.XCLOSE)THEN
        ILOC=I
        VLOC=TMP(I)
        XCLOSE=X
       ENDIF 
20    CONTINUE
C
      RETURN
      END
