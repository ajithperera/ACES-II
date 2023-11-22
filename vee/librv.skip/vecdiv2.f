      SUBROUTINE VECDIV2(ROOT,A,B,C,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(N),C(N)
      DO 10 I=1,N
       X=ROOT-B(I)
       IF(ABS(X).LT.1.D-4)THEN
        X=SIGN(1.0D-4,X)
       ENDIF
       C(I)=A(I)/X
10    CONTINUE
      RETURN
      END
