      SUBROUTINE MHOUSE(NN,MM,MM1,G,X,A)
C  ....... This way didn't work ......................
C     SUBROUTINE MHOUSE(NN,MM,G,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8           SIGMA,SUM,AI4,GIJ,XIJ,ANORM,MU,LAMBDA
C     DIMENSION G(MM,1),X(MM,1),A(MM,5)
C  ........... I am not sure about the size of A matrix
C     DIMENSION G(MM,1),X(MM,1),A(300,5)
      DIMENSION G(MM,1),X(MM,1),A(MM1,5)
C  ....... MM1 = MM + 1, MM.GE.NN
C   TOLER IS THE RELATIVE PRECISION OF THE MACHINE
      TOLER=1.D-15
      N=NN
      N1=N-1
C  ....  This is a kind of stupid limitation of dimension
C         IF(N.GT.MM .OR.N.GT.100)  GO TO 21
      IF(N1) 21,24,23
   24 X(1,1)=1.D0
      GO TO 21
C SET XTO IDENTITY
   23 DO 1 I=1,N
      X(I,I)=1.D0
      IF(I.EQ.N) GO TO 22
      I1=I+1
      DO 1 J=I1,N
      X(I,J)=0.D0
    1 X(J,I)=0.D0
C  REDUCE G TO TRIDIAGONAL FORM USING N-2 OTHOGONAL TRANSFORMS
C       (I-2WW) = ( I-GAMMA*UU)
C  G IS REPLACED BY U BELOW THE MAIN DIAGONAL.
C  NORM IS THE INFINITY NORM OF THE REDUCED MATRIX.
   22 ANORM=0.0D0
      ABSB=0.0D0
      N2=N-2
      IF(N2.LE.0) GO TO 25
      DO 2 K=1,N2
      A(K,1)=G(K,K)
      K1=K+1
      SIGMA=0.0D0
      DO 3 I=K1,N
      GG=G(I,K)**2
    3 SIGMA=SIGMA+GG
      T=ABSB+DABS(A(K,1))
      ABSB=DSQRT(SIGMA)
      TA=T+ABSB
      ANORM=DMAX1(ANORM,TA)
      ALPHA=G(K1,K)
      BETA=ABSB
      IF(ALPHA.GE.0.0) BETA=-BETA
      A(K1,2)=BETA
      IF(SIGMA.EQ.0.0D0) GO TO 2
      GAMMA=1.0D0/(SIGMA-ALPHA*BETA)
      G(K1,K)=ALPHA-BETA
      DO 6 I=K1,N
      SUM=0.0D0
      DO 7 J=K1,I
       GG=G(I,J)*G(J,K)
    7 SUM=SUM+GG
      IF(I.EQ.N) GO TO 6
      I1=I+1
      DO 8 J=I1,N
      GG=G(J,I)*G(J,K)
    8 SUM=SUM+GG
    6 A(I,4)=SUM*GAMMA
      SUM=0.0D0
      DO 9 I=K1,N
      GG=G(I,K)*A(I,4)
    9 SUM=SUM+GG
      T=0.5D0*GAMMA*SUM
      DO 10 I=K1,N
      AI4=A(I,4)
   10 A(I,4)=AI4-T*G(I,K)
      DO 11 I=K1,N
      DO 11 J=K1,I
      GIJ=G(I,J)
   11 G(I,J)=GIJ-G(I,K)*A(J,4)-A(I,4)*G(J,K)
      DO 12 I=2,N
      SUM=0.0D0
      DO 13 J=K1,N
      GA=X(I,J)
      GB=G(J,K)
   13 SUM=SUM+GA*GB
   12 A(I,4)=GAMMA*SUM
      DO 14 I=2,N
      DO 14 J=K1,N
      XIJ=X(I,J)
   14 X(I,J)=XIJ-A(I,4)*G(J,K)
    2 CONTINUE
   25 CONTINUE
      A(N1,1)=G(N1,N1)
      A(N,1)=G(N,N)
      A(N,2)=G(N,N1)
      T=DABS(A(N,2))
      TA=DABS(A(N1,1))+T+ABSB
      TAA= T+DABS(A(N,1))
      ANORM=DMAX1(ANORM,TA,TAA)
      EPS=ANORM*TOLER
C  INFINITY NORM TIMES RELATIVE MACHINE PRECISION
      MU=0.D0
      A(1,2)=0.D0
      A(N+1,2)=0.D0
      M=N
   15 IF(M.EQ.0) GO TO 21
      I=M-1
      M1=I
      K=I
      IF(DABS(A(K+1,2)).GT.EPS) GO TO 16
      G(M,M)=A(M,1)
      MU=0.D0
      M=K
      GO TO 15
   16 IF(DABS(A(I,2)).LE.EPS.OR.I.EQ.1) GO TO 18
      I=I-1
      K=I
      GO TO 16
   18 LAMBDA=0.D0
      IF(DABS(A(M,1)-MU).LT.0.5D0*DABS(A(M,1)).OR.M1.EQ.K)
     1    LAMBDA=A(M,1)+0.5D0*A(M1+1,2)
      MU=A(M,1)
      A(K,1)=A(K,1)-LAMBDA
      BETA=A(K+1,2)
C  TRANSFORMATION ON THE LEFT
      DO 19 J=K,M1
      A0=A(J,1)
      A1=A(J+1,1)-LAMBDA
       B0=A(J+1,2)
      T=DSQRT(A0**2+BETA**2)
      COSINE=A0/T
      A(J,4)=COSINE
      SINE=BETA/T
      A(J,5)=SINE
      A(J,1)=COSINE*A0+SINE*BETA
      A(J+1,1)=-SINE*B0+COSINE*A1
      A(J+1,2)=COSINE*B0+SINE*A1
      BETA=A(J+2,2)
      A(J+2,2)= COSINE*BETA
   19 A(J+2,3)=SINE*BETA
      A(K,2)=0.0
      A(K+1,3)=0.0
C  TRANSFORMATION ON THE RIGHT
      DO 20 J=K,M1
      SINE=A(J,5)
      COSINE=A(J,4)
      A0=A(J,1)
      B0=A(J+1,2)
      A(J,2)=A(J,2)*COSINE+A(J+1,3)*SINE
      A(J,1)=A0*COSINE+B0*SINE +LAMBDA
      A(J+1,2)=-A0*SINE+B0*COSINE
      A(J+1,1)=A(J+1,1)*COSINE
      DO 20 I=1,N
      X0=X(I,J)
      X1=X(I,J+1)
      X(I,J)=X0*COSINE+X1*SINE
   20 X(I,J+1)=-X0*SINE+X1*COSINE
      A(M,1)=A(M,1)+LAMBDA
      GO TO 15
   21 CONTINUE
      RETURN
      END
