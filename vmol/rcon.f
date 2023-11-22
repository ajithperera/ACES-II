      SUBROUTINE RCON(H2,BX,BY,BZ,H,N,J1,J2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BX(N,J1),BY(N,J2),BZ(N,J1+J2),H(N,2),H2(2)
      IF(N .GT. MAX0(J1,J2)) THEN
      DO 100 JJ=1,J2
      DO 1 K=1,N
    1 H(K,JJ) = BX(K,1)*BZ(K,JJ+1)
  100 CONTINUE
      DO 2 II=2,J1
      DO 2 JJ=1,J2
      DO 2 K=1,N
    2 H(K,JJ) = H(K,JJ) + BX(K,II)*BZ(K,II+JJ)
      DO 3 K=1,N
    3 H2(K) = H(K,1)*BY(K,1)
      DO 4 JJ=2,J2
      DO 4 K=1,N
    4 H2(K) = H2(K) + H(K,JJ)*BY(K,JJ)
C
      ELSE IF (J1 .LT. J2) THEN
      DO 5 K=1,N
      DO 5 JJ=1,J2
    5 H(K,JJ) = BX(K,1)*BZ(K,JJ+1)
      DO 6 II=2,J1
      DO 6 K=1,N
      DO 6 JJ=1,J2
    6 H(K,JJ) = H(K,JJ) + BX(K,II)*BZ(K,II+JJ)
      DO 7 K=1,N
    7 H2(K) = H(K,1)*BY(K,1)
      DO 8 JJ=2,J2
      DO 8 K=1,N
    8 H2(K) = H2(K) + H(K,JJ)*BY(K,JJ)
C
      ELSE
      DO 9 K=1,N
      DO 9 JJ=1,J1
    9 H(K,JJ) = BY(K,1)*BZ(K,JJ+1)
      DO 10 II=2,J2
      DO 10 K=1,N
      DO 10 JJ=1,J1
   10 H(K,JJ) = H(K,JJ) + BY(K,II)*BZ(K,II+JJ)
      DO 11 K=1,N
   11 H2(K) = H(K,1)*BX(K,1)
      DO 12 JJ=2,J1
      DO 12 K=1,N
   12 H2(K) = H2(K) + H(K,JJ)*BX(K,JJ)
      ENDIF
      RETURN
      END
