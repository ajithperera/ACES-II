      SUBROUTINE GAUSSL(X1,X2,N,X,W)
C
C THIS ROUTINE DETERMINES THE ROOTS AND WEIGHTS
C FOR THE GAUSS-LEGENDRE INTEGRATION.       
C
C INPUT: X1 ..... LOWER LIMIT OF INTEGRATION RANGE
C        X2 ..... UPPER LIMIT OF INTEGRATION RANGE
C        N  ..... NUMBER OF INTEGRATION POINTS
C
C OUTPUT X  ..... ROOTS
C        W  ..... WEIGHTS
C
CEND
C
C JG 4/93, TAKEN FROM NUMERICAL RECIPE.
C      
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),W(N)
      PARAMETER (EPS=3.D-14)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
       Z=COS(3.141592654D0*(I-.25D0)/(N+0.5D0))
1      CONTINUE
        P1=1.D0
        P2=0.D0
        DO 11 J=1,N
         P3=P2
         P2=P1
         P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11      CONTINUE
        PP=N*(Z*P1-P2)/(Z*Z-1.D0)
        Z1=Z
        Z=Z1-P1/PP
       IF(ABS(Z-Z1).GT.EPS) GO TO 1
       X(I)=XM-XL*Z
       X(N+1-I)=XM+XL*Z
       W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
       W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END