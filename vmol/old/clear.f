      SUBROUTINE CLEAR (A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N)
      DO 1 I=1,N
    1 A(I)=0.
      RETURN
      END
