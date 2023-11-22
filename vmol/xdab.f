
      SUBROUTINE XDAB(X,D,A,B,N)
      IMPLICIT NONE
      DOUBLE PRECISION X,D,A,B
      INTEGER N,I
      DIMENSION X(N),A(N),B(N)
C
C     This routine does a type of vector product :
C
C     X(I) = X(I) + D*A(I)*B(I)
C
C     If there is something that does the same thing in blas, we'll use
C     that, but for now this will have to suffice.
C
      DO 10 I=1,N
      X(I) = X(I) + D * A(I) * B(I)
   10 CONTINUE
      RETURN
      END
