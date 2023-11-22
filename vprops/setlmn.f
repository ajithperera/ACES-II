C
      SUBROUTINE SETLMN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXANG=7)
C
C....    SET ANGULAR QUANTUM NUMBERS AND COMPUTE CARTESIAN
C....    NORMALIZATION FACTORS
C
      COMMON /HIGHL/ LMNVAL(3,84), ANORM(84)
      DIMENSION FACS(0:MAXANG-1)
      DATA FACS /1.D0, 1.D0, 3.D0, 15.D0, 105.D0, 945.D0, 10395.D0/
      II = 0
      DO 10 LVAL = 0,MAXANG-1
         DO 20 L = LVAL,0,-1
            LEFT = LVAL - L
            DO 30 M = LEFT,0,-1
               N = LEFT - M
               II = II + 1
               LMNVAL(1,II) = L
               LMNVAL(2,II) = M
               LMNVAL(3,II) = N
               ANORM(II) = FACS(L)*FACS(M)*FACS(N)
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
      RETURN 
      END
