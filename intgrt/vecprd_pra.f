
C COMPUTES C(I)=A(I)*B(I) FOR FIRST N ELEMENTS OF VECTORS A AND B.

      SUBROUTINE VECPRD_PRA(A,B,C,N,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,M),B(N,M),C(N,M)
      if (n.lt.1) return
      DO I = 1, N
        do j=1,m
         C(I,j) = A(I,j) * B(I,j)
        end do
      END DO
      RETURN
      END
