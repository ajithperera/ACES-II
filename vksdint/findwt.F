      SUBROUTINE findwt(X1,X2,X,W,N)
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (EPS =3.0D-14)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE =1.0D+00)
      PARAMETER (TWO =2.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (PT5 =0.5D+00)
      DIMENSION X(*),W(*)
      PI=FOUR*ATAN(ONE)
      NR=(N+1)/2
      XM=PT5*(X2+X1)
      XL=PT5*(X2-X1)
      PIN=(PI/TWO)/(2*N+1)
      DO I=1,NR
         Z=COS(PIN*(4*I-1))
    1    CONTINUE
            P1=ONE
            P2=ZERO
            DO J=1,N
               P3=P2
               P2=P1
               P1=((J+J-1)*Z*P2-(J-1)*P3)/J
            ENDDO
            PP=N*(Z*P1-P2)/(Z*Z-ONE)
            Z1=Z
            Z =Z1-P1/PP
            IF(ABS(Z-Z1).GT.EPS) GO TO 1
         X(    I)=XM+XL*Z
         X(N+1-I)=XM-XL*Z
         W(    I)=TWO*XL/((ONE-Z*Z)*PP*PP)
         W(N+1-I)=W(I)
      ENDDO
      RETURN
      END

