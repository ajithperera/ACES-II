      SUBROUTINE FNDCLOSE(LENGTH,V,TEST,VLOCMIN,ILOCMIN)
C
C LOCATES THE ELEMENT IN A VECTOR WHICH IS CLOSEST TO TEST 
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(LENGTH)
C
      XCLOSE=1.D+30
      DO 10 I=1,LENGTH
       TMP=ABS(V(I)-TEST)
       IF(TMP.LT.XCLOSE)THEN
        ILOCMIN=I
        VLOCMIN=V(I)
        XCLOSE=TMP
       ENDIF 
10    CONTINUE
C
      RETURN
      END
