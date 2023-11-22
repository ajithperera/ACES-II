      SUBROUTINE HGSOL (NDIM, N, ML, MU, B, Y, IP)
C
C
C     Description of routine.
C
C     Solution of  A*X = C  given LU decomposition of a from HGDEC.
C
C     Y  =  right-hand vector C, of length N, on input,
C
C        =  solution vector X on output.
C
C     All the arguments are input arguments.
C     The output argument is  Y.
C
C     Original code: copyright A. C. Hindmarsh and LLL, 1975.
C     Modifications: copyright Erik Deumens and QTP, 1986, 1992.
C
      IMPLICIT NONE

      INTEGER NDIM, N
      DOUBLE PRECISION
     $   B(NDIM,1),Y(N), XX, DP
      INTEGER IP(N)
      INTEGER ML, MU, N1, LL, NR, J, KK, I, NB
      IF (N .EQ. 1) GO TO 60
      N1 = N - 1
      LL = ML + MU + 1
      IF (ML .EQ. 0) GO TO 32
      DO 30 NR = 1,N1
        IF (IP(NR) .EQ. NR) GO TO 10
        J = IP(NR)
        XX = Y(NR)
        Y(NR) = Y(J)
        Y(J) = XX
 10     KK = MIN0(N-NR,ML)
        DO 20 I = 1,KK
 20       Y(NR+I) = Y(NR+I) + Y(NR)*B(NR,LL+I)
 30     CONTINUE
 32   LL = LL - 1
      Y(N) = Y(N)*B(N,1)
      KK = 0
      DO 50 NB = 1,N1
        NR = N - NB
        IF (KK .NE. LL) KK = KK + 1
        DP = 0.D0
        IF (LL .EQ. 0) GO TO 50
        DO 40 I = 1,KK
 40       DP = DP + B(NR,I+1)*Y(NR+I)
 50     Y(NR) = (Y(NR) - DP)*B(NR,1)
      RETURN
 60   Y(1) = Y(1)*B(1,1)
      RETURN
      END
