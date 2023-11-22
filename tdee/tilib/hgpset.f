      SUBROUTINE HGPSET (DIFFUN, JAC, ISTATS, Y, N0, CON, MITER, IER,
     $      YMAX, ERROR, SAVE1, SAVE2, PW, IPIV)
C
C
C     Description of routine.
C     HGPSET is called by HGSTEP to compute and process the matrix
C     P = I - H*EL(1)*J , where J is an approximation to the jacobian.
C     This is the version for banded form of J.
C     J is computed, either by the user-supplied routine JAC
C     if MITER = 1, or by finite differencing if MITER = 2.
C     J is stored in PW and replaced by P, using CON = -H*EL(1).
C     Then P is subjected to LU decomposition in preparation for
C     later solution of linear systems with P as coefficient matrix.
C
C     In addition to variables described previously, communication
C     with HGPSET uses the following:
C
C     EPSJ    = SQRT(UROUND), used in the numerical jacobian increments.
C
C     MW      = ML + MU + 1.
C
C     NM1     = N0 - 1.
C
C     N0ML    = N0*ML.
C
C     N0W     = N0*MW.
C
C     Original code: copyright A. C. Hindmarsh and LLL, 1975.
C     Modifications: copyright Erik Deumens and QTP, 1986, 1992.
C
      IMPLICIT NONE

      INTEGER N0
      DOUBLE PRECISION Y(N0,1), CON
      EXTERNAL DIFFUN, JAC
      DOUBLE PRECISION T, H, HMIN,HMAX,EPS,UROUND
      INTEGER N, MFC,KFLAG,JSTART, MITER, IER, ISTATS
      COMMON /GEAR1/ T,H,HMIN,HMAX,EPS,UROUND,N,MFC,KFLAG,JSTART
C     Note:      ED QTP 31-Dec-1986
C     Arrays are passed as arguments to make them dynamic
C      COMMON /GEAR2/ YMAX(1)
C      COMMON /GEAR3/ ERROR(1)
C      COMMON /GEAR4/ SAVE1(1)
C      COMMON /GEAR5/ SAVE2(1)
C      COMMON /GEAR6/ PW(1)
C      COMMON /GEAR7/ IPIV(1)
      DOUBLE PRECISION
     $   YMAX(1), ERROR(1), SAVE1(1), SAVE2(1), PW(1), EPSJ,
     $   R0, D, R
      INTEGER IPIV(1), I, II, I1, I2, J, J1, JJ, K, KMAX
      INTEGER ML, MU, MW, NM1, N0ML, N0W
      COMMON /GEAR8/ EPSJ,ML,MU,MW,NM1,N0ML,N0W
C-----------------------------------------------------------------------
      IF (MITER .EQ. 2) GO TO 20
C     If MITER = 1, call JAC and multiply by scalar.
      CALL JAC (N, T, Y, PW, N0, ML, MU, ISTATS)
      IF (MOD(ISTATS,2) .NE. 1) RETURN
      DO 10 I = 1,N0W
 10     PW(I) = PW(I)*CON
      GO TO 90
C     If MITER = 2, make MW calls to DIFFUN to approximate J.
 20   D = 0.D0
      DO 30 I = 1,N
 30     D = D + SAVE2(I)**2
      R0 = ABS(H)*SQRT(D)*1.D03*UROUND
C     The original values of Y(,1) are saved temporarily in ERROR().
      DO 40 I = 1,N
 40     ERROR(I) = Y(I,1)
      DO 80 J = 1,MW
        KMAX = (N - J)/MW + 1
        DO 50 K = 1,KMAX
          I = J + (K - 1)*MW
          R = MAX(EPSJ*YMAX(I),R0)
 50       Y(I,1) = Y(I,1) + R
        CALL DIFFUN (T, Y, SAVE1, ISTATS)
        IF (MOD(ISTATS,2) .NE. 1) RETURN
        J1 = J*N0 + N0ML
        DO 70 K = 1,KMAX
          JJ = J + (K - 1)*MW
          R = MAX(EPSJ*YMAX(JJ),R0)
          D = CON/R
          I1 = MAX0(JJ-MU,1)
          I2 = MIN0(JJ+ML,N)
          II = J1 - NM1*I1
          DO 60 I = I1,I2
            PW(II) = (SAVE1(I) - SAVE2(I))*D
 60         II = II - NM1
          Y(JJ,1) = ERROR(JJ)
          ERROR(JJ) = 0.D0
 70       J1 = J1 + N0W
 80     CONTINUE
C     Add identity matrix.
 90   DO 100 I = 1,N
 100    PW(N0ML+I) = PW(N0ML+I) + 1.D0
C     Do LU decomposition on P.
      CALL HGDEC (N0, N, ML, MU, PW, IPIV, IER)
      RETURN
      END
