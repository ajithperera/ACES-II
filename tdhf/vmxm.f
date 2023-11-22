C                                                                  ************
      SUBROUTINE VMXM(A,B,C,N)
C
C
      INTEGER N,I,J,K
      DOUBLE PRECISION A(N,N),B(N),C(N)
C
         DO 20 I = 1,N
            C(I) = 0.0
            DO 10 K = 1,N
               C(I) = C(I) + A(I,K)*B(K)
   10       CONTINUE
   20    CONTINUE
      RETURN
      END
