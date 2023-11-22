            SUBROUTINE EIG(A,B,L,N,N1)
C
C DIAGONALIZATION BY THE JACOBI METHOD.
C A - MATRIX TO BE DIAGONALIZED
C B - EIGENVECTORS
C L - DIMENSION OF A AND B
C N - SIZE OF SUBMATRIX USED
C N1 - A FLAG INDICATING WHETHER THE EIGENVECTORS AND
C      EIGENVALUES ARE TO BE REORDERED.
C 
CSW      1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(L,L),B(L,L)
      DATA ZERO/0.D00/,ONE/1.D00/
CSW     2
      ABS(X)=DABS(X)
      SQRT(X)=DSQRT(X)
      TOL=1.D-14
      B(1,1)=ONE
      IF(N.EQ.1) RETURN
      DO 1 I=2,N
      IM1=I-1
      DO 2 J=1,IM1
      B(I,J)=ZERO
    2 B(J,I)=ZERO
    1 B(I,I)=ONE
  200 P=ZERO
      DO 3 I=2,N
      IM1=I-1
      DO 4 J=1,IM1
      Q=A(I,J)
      IF(P.GE. ABS(Q)) GO TO 4
      P= ABS(Q)
      II=I
      JJ=J
    4 CONTINUE
    3 CONTINUE
      IF(P.EQ.0.) GO TO 205
      P=A(II,II)
      Q=A(II,JJ)
      R=A(JJ,JJ)
      DIFF=0.5D 00*(P-R)
      IF( ABS(DIFF).LT. ABS(Q)) GO TO 201
      IF( ABS(Q/DIFF).GT.TOL) GO TO 201
      A(II,JJ)=ZERO
      A(JJ,II)=ZERO
      GO TO 200
  201 S=SQRT(0.250*(P-R)**2+Q**2)
      SUM=0.5D 00*(P+R)
      D=R*P-Q**2
      IF(SUM.GT.ZERO) GO TO 5
      ALN=SUM-S
      ALP=D/ALN
      GO TO 6
    5 ALP=SUM+S
      ALN=D/ALP
    6 IF(DIFF.GT.ZERO) GO TO 7
      T=Q/(DIFF-S)
      A(II,II)=ALN
      A(JJ,JJ)=ALP
      GO TO 8
    7 T=Q/(DIFF+S)
      A(II,II)=ALP
      A(JJ,JJ)=ALN
    8 C=1.0/SQRT(1.0+T**2)
      S=T*C
      A(II,JJ)=ZERO
      A(JJ,II)=ZERO
      DO 9 I=1,N
      P=B(I,II)
      Q=B(I,JJ)
      B(I,II)=C*P+S*Q
      B(I,JJ)=C*Q-S*P
      IF(I.EQ.II.OR.I.EQ.JJ) GO TO 9
      P=A(I,II)
      Q=A(I,JJ)
      R=C*P+S*Q
      A(I,II)=R
      A(II,I)=R
      R=Q*C-P*S
      A(I,JJ)=R
      A(JJ,I)=R
    9 CONTINUE
      GO TO 200
  205 IF(N1.EQ.1) RETURN
      MM=N-1
      DO 20 I=1,MM
      II=I+1
      DO 20 J=II,N
      IF(A(I,I)-A(J,J)) 15,15,20
   15 W1=A(I,I)
      A(I,I)=A(J,J)
      A(J,J)=W1
      DO 10 MU=1,N
      W2=B(MU,I)
      B(MU,I)=B(MU,J)
   10 B(MU,J)=W2
   20 CONTINUE
      RETURN
      END
