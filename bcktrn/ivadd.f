 
      SUBROUTINE IVADD(A,B,C,N)
      IMPLICIT INTEGER (A-Z)
      DIMENSION A(N),B(N),C(N)
      DO 10 I=1,N
       C(I)=A(I)+B(I)
10    CONTINUE 
      RETURN
      END
