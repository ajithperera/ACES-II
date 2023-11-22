      SUBROUTINE COLLEC (K,A,IB,M,TMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M),IB(M)
      K=1
      DO 21 I=1,M
      IB(K)=I
c
cjfs -original code.  good machine independent stuff.
c      K=K+ISHFT(A(I)-TMAX,63)
c
c replacement code
c
       Z=A(I)-TMAX
       IF(Z.LT.0.0)K=K+1
   21 CONTINUE
      K=K-1
      RETURN
      END
