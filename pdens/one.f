
C ZERO OUT THE FIRST LEN ELEMENTS OF DOUBLE PRECISION VECTOR A.

      SUBROUTINE ONE(A,LEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(LEN)
      if (len.lt.1) return
      DO I = 1, LEN
         A(I) = 1.0D0*A(I)
      END DO
      RETURN
      END
