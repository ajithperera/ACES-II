      SUBROUTINE VECDIV2(ROOT,A,B,C,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PRINT
      DIMENSION A(N),B(N),C(N)
      DATA ZILCH /0.0D0/
      print = .false.
      DO 10 I=1,N
       X=ROOT-B(I)
       IF (print .and. (ABS(X).lt.0.005)) then
          write(6,*)'  small denominator ', x, i
       endif
       IF(ABS(X).LT.1.D-4)THEN
        X=SIGN(1.0D-4,X)
       ENDIF
       C(I)=A(I)/X
10    CONTINUE
      RETURN
      END
