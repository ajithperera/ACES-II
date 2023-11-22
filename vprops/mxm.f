      SUBROUTINE MXM(A,LL,B,MM,C,NN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER L, M, N
      DIMENSION A(LL,MM), B(MM,NN), C(LL,NN)
      DATA ONE /1.0D0/
      DATA ZILCH /0.0D0/
      L = LL
      M = MM
      N = NN
      CALL XGEMM('N','N',L,N,M,ONE,A,L,B,M,ZILCH,C,L)
      RETURN
      END
