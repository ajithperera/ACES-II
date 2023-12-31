
      SUBROUTINE FNDMINE(LENGTH,VINPUT,VEXCLUDE,TMP,NEXCLUDE,
     &                   VLOCMIN,ILOCMIN,TOL)
C
C FINDS MINIMUM ELEMENT IN VECTOR VINPUT WHICH IS NOT WITHIN TOL
C OF THE VALUES VEXCLUDE(1..NEXCLUDE).
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VINPUT(LENGTH),VEXCLUDE(NEXCLUDE),TMP(LENGTH)
C
      CALL SCOPY(LENGTH,VINPUT,1,TMP,1)
      DO 10 I=1,NEXCLUDE
       DO 11 J=1,LENGTH
        CALL FNDCLOSE(LENGTH,TMP,VEXCLUDE(I),X,ILOC)
        X2=ABS(VEXCLUDE(I)-TMP(ILOC))
        IF(X2.LT.TOL)TMP(ILOC)=1.D+30
11     CONTINUE
10    CONTINUE
C
      VLOCMIN=1.D+30
      DO 15 I=1,LENGTH
       X=TMP(I)
       IF(X.LT.VLOCMIN)THEN
        ILOCMIN=I
        VLOCMIN=TMP(I)
       ENDIF
15    CONTINUE
C
      RETURN
      END
