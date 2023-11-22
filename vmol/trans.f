      SUBROUTINE TRANS (NLOP,C,ILX,X,NN,MM,ISTO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(2),ILX(2),X(2),NN(2),MM(2)
C
      MMM=1
      DO 1 L=1,NLOP
    1 MMM=MMM*NN(L)
      ISTO=0
C     LOOP OVER INDICES TO BE TRANSFORMED
      DO 6 L=1,NLOP
      N=NN(L)
      M=MM(L)
      IL = ILX(L)
      MN=MMM/N
        call mxm(x,mn,c(il+1),n,x(mmm+1),m)
c        call trapo(x,x(mmm+1),m,mn)
        call transp(x(mmm+1),x,m,mn)
    6 MMM=MN*M
      RETURN
      END
