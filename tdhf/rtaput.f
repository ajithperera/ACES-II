      SUBROUTINE RTAPUT(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N),B(1)
      IJ=0
      DO 1 I=1,N
      DO 1 J=1,I
      IJ=IJ+1
      A(I,J)=B(IJ)
      A(J,I)=-A(I,J)
    1 CONTINUE
      RETURN
      END
