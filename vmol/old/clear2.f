      SUBROUTINE CLEAR2 (A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),B(N)
      DO 2 I=1,N
      A(I)=0.
    2 B(I)=0.
      RETURN
      END
