      SUBROUTINE IDXLUK(N,X)
C
      IMPLICIT INTEGER(A-Z)
C
      DIMENSION X(N*(N+1)/2,2)
C
      DO 100 I=1,N
        DO 110 J=1,I
          INDEX=J+I*(I-1)/2
          X(INDEX,1)=I
          X(INDEX,2)=J
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
