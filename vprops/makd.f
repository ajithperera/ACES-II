      SUBROUTINE MAKD(DMAT,VEC,NAO,NMO,FAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DMAT(1),VEC(NAO,NMO),FAC(1)
      IX=0
      DO 1 I=1,NAO
      DO 2 J=1,I
      X=0.0
      DO 3 K=1,NMO
      X=X+VEC(I,K)*VEC(J,K)*FAC(K)
3     CONTINUE
      IX=IX+1
      IF(I.NE.J)X=X*2.0
      DMAT(IX)=X
2     CONTINUE
1     CONTINUE
      RETURN
      END