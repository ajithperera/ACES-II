      SUBROUTINE MXMT(A,LL,B,MM,C,NN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C....   DOES C = AT*B
C
      INTEGER L, M, N
      DIMENSION A(LL,MM), B(MM,NN), C(LL,NN)
      DATA ONE /1.0D0/
      DATA ZILCH /0.0D0/
      L = LL
      M = MM
      N = NN
      CALL XGEMM('T','N',L,N,M,ONE,A,M,B,M,ZILCH,C,L)
      RETURN
      END
