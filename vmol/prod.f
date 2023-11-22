      SUBROUTINE PROD(A,B,C,M,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N) , B(M) , C(M,N)
      IF (M .GE. N)THEN
      DO 1 I=1,N
      DO 1 J=1,M
    1 A(J,I) = B(J) * C(J,I)
C
      ELSE
      DO 2 J=1,M
      DO 2 I=1,N
    2 A(J,I) = B(J) * C(J,I)
      ENDIF
      RETURN
      END
