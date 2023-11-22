      SUBROUTINE TRANSQ(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N)
      DO 10 I=1,N
      DO 10 J=1,I
      X=A(I,J)
      A(I,J)=A(J,I)
      A(J,I)=X
   10 CONTINUE
      RETURN
      END 
