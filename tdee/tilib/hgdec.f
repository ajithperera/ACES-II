      SUBROUTINE HGDEC (NDIM, N, ML, MU, B, IP, IER)
C
C
C     Description of routine.
C
C     LU decomposition of band matrix A: L*U = P*A , where P is a
C        permutation matrix, L is a unit lower triangular matrix,
C        and U is an upper triangular matrix.
C
C     N     =  Order of matrix.
C
C     B     =  N by (2*ML+MU+1) array containing the matrix A on input
C              and its factored form on output.
C              On input, B(I,K) (1.LE.I.LE.N) contains the K-th
C              diagonal of A, or A(I,J) is stored in B(I,J-I+ML+1).
C              On output, B contains the L and U factors, with
C              U in columns 1 to ML+MU+1, and L in columns
C              ML+MU+2 to 2*ML+MU+1.
C
C     ML,MU =  Widths of the lower and upper parts of the band, not
C              counting the main diagonal. Total bandwidth is ML+MU+1.
C
C     NDIM  =  The first dimension (column length) of the array B.
C              NDIM must be .GE. N.
C
C     IP    =  Array of length N containing pivot information.
C
C     IER   =  Error indicator:
C
C           =  0  if no error,
C
C           =  K  if the K-th pivot chosen was zero (A is singular).
C
C     The input arguments are  NDIM, N, ML, MU, B.
C
C     The output arguments are  B, IP, IER.
C
C     Original code: copyright A. C. Hindmarsh and LLL, 1975.
C     Modifications: copyright Erik Deumens and QTP, 1986, 1992.
C
      IMPLICIT NONE

      INTEGER ML, MU, NDIM, IER, N, I, K, II, NR, N1, MX, NP
      DOUBLE PRECISION B(NDIM,1), XM, XX
      INTEGER IP(N)
      INTEGER LL, J, KK, LR
C-----------------------------------------------------------------------
      IER = 0
      IF (N .EQ. 1) GO TO 92
      LL = ML + MU + 1
      N1 = N - 1
      IF (ML .EQ. 0) GO TO 32
      DO 30 I = 1,ML
        II = MU + I
        K = ML + 1 - I
        DO 10 J = 1,II
 10       B(I,J) = B(I,J+K)
        K = II + 1
        DO 20 J = K,LL
 20       B(I,J) = 0.D0
 30     CONTINUE
 32   LR = ML
      DO 90 NR = 1,N1
        NP = NR + 1
        IF (LR .NE. N) LR = LR + 1
        MX = NR
        XM = ABS(B(NR,1))
        IF (ML .EQ. 0) GO TO 42
        DO 40 I = NP,LR
          IF (ABS(B(I,1)) .LE. XM) GO TO 40
          MX = I
          XM = ABS(B(I,1))
 40       CONTINUE
 42     IP(NR) = MX
        IF (MX .EQ. NR) GO TO 60
        DO 50 I = 1,LL
          XX = B(NR,I)
          B(NR,I) = B(MX,I)
 50       B(MX,I) = XX
 60     XM = B(NR,1)
        IF (XM .EQ. 0.D0) GO TO 100
        B(NR,1) = 1.D0/XM
        IF (ML .EQ. 0) GO TO 90
        XM = -B(NR,1)
        KK = MIN0(N-NR,LL-1)
        DO 80 I = NP,LR
          J = LL + I - NR
          XX = B(I,1)*XM
          B(NR,J) = XX
          DO 70 II = 1,KK
 70         B(I,II) = B(I,II+1) + XX*B(NR,II+1)
 80       B(I,LL) = 0.D0
 90     CONTINUE
 92   NR = N
      IF (B(N,1) .EQ. 0.D0) GO TO 100
      B(N,1) = 1.D0/B(N,1)
      RETURN
 100  IER = NR
      RETURN
      END
