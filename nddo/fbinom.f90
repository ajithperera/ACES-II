      SUBROUTINE FBINOM
      use tables
!     *
!     DEFINE FACTORIALS AND BINOMIAL COEFFICIENTS (FOR SPD INTEGRALS).
!     F(30)         LOGARITHM OF FACTORIALS.
!     B(30,30)      BINOMIAL COEFFICIENTS.
!     *
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      COMMON
!     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
!     ./PSC   / F(30),B(30,30)
      DATA ZERO  /  0.0D0/
      DATA ONE   /  1.0D0/
      DATA TWO   /  2.0D0/
      DATA THREE /  3.0D0/
      DATA FOUR  /  4.0D0/
      DATA PT5   /  0.5D0/
      DATA PT25  /  0.25D0/
      K      = 30
      F(1)   = ZERO
      DO 20 I=1,29
   20 F(I+1) = F(I)+LOG(ZERO+I)
      DO 30 I=1,K
      B(I,1) = ONE
      DO 30 J=2,K
   30 B(I,J) = ZERO
      DO 40 I=2,K
      DO 40 J=2,I
   40 B(I,J) = B(I-1,J-1)+B(I-1,J)
      RETURN
      END subroutine fbinom
