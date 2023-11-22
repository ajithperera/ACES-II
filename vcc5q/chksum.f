
C SUMS THE FIRST LEN ELEMENTS OF VECTOR V AND WRITES IT TO STANDARD OUTPUT.

      SUBROUTINE CHKSUM(V,LEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision V(*)
      Z = 0.0d0
      DO I = 1, LEN
         Z = Z + V(I)
      END DO
      WRITE(*,*) Z
      RETURN
      END
