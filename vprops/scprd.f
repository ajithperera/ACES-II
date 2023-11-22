      FUNCTION SCPRD(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),B(1)
      X=0.0
      DO 1 I=1,N
      X=X+A(I)*B(I)
1     CONTINUE
      SCPRD=X
      RETURN
      END
