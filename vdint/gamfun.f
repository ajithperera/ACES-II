
      SUBROUTINE GAMFUN
C
C     Trygve Ulf Helgaker fall 1984
C
C     This subroutine calculates the incomplete gamma function as
C     described by McMurchie & Davidson, J. Comp. Phys. 26 (1978) 218.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (ONE = 1.D0,  TWO = 2.D0, TEN = 10.D0,
     *           HALF = .5D0, TENTH = .1D0)
      PARAMETER (SQRPIH = 0.88622 69254 52758 D00)
      PARAMETER (COEF2 = 0.5000000000000000 D00,
     *           COEF3 = -.1666666666666667 D00,
     *           COEF4 = 0.0416666666666667 D00,
     *           COEF5 = -.0083333333333333 D00,
     *           COEF6 = 0.0013888888888889 D00)
      PARAMETER (GFAC30 =  .4999489092 D0,
     *           GFAC31 = -.2473631686 D0,
     *           GFAC32 =  .321180909  D0,
     *           GFAC33 = -.3811559346 D0,
     *           GFAC20 =  .4998436875 D0,
     *           GFAC21 = -.24249438   D0,
     *           GFAC22 =  .24642845   D0,
     *           GFAC10 =  .499093162  D0,
     *           GFAC11 = -.2152832    D0,
     *           GFAC00 =  .490        D0)
C
C     GAMCOM: 3267 = 27*121, max J value = 20
      PARAMETER (MXQN=8)
      PARAMETER (MXAQN=MXQN*(MXQN+1)/2,MXAQNS=MXAQN*MXAQN)
      PARAMETER (MAXJ = 32)
CSSS      COMMON /GAMCOM/ WVAL, FJW(0:MAXJ), TABFJW(3267), JMAX0
      COMMON/GAMCOM/WVAL,FJW(0:4*(MXQN-1)+2), 
     *              TABFJW(121*(4*(MXQN-1)+2+7)), JMAX0
C
      SAVE MAXJ0
      DATA MAXJ0 /-1/
C
      IF (JMAX0 .GT. MAXJ0) THEN
         WRITE (LUPRI,'(//A,I5,A,I3,A)')
     *     ' GAMFUN ERROR: JMAX0 =',JMAX0,', which is greater than',
     *     MAXJ0,' from GAMTAB call.'
         STOP 
      END IF
      IPOINT = IDNINT(TEN*WVAL)
      IF (IPOINT .LT. 120) THEN
         ISTART = 1 + 121*JMAX0 + IPOINT
         WDIF = WVAL - TENTH*FLOAT(IPOINT)
         FJW(JMAX0) = (((((COEF6*TABFJW(ISTART + 726)*WDIF
     *                   + COEF5*TABFJW(ISTART + 605))*WDIF
     *                    + COEF4*TABFJW(ISTART + 484))*WDIF
     *                     + COEF3*TABFJW(ISTART + 363))*WDIF
     *                      + COEF2*TABFJW(ISTART + 242))*WDIF
     *                       - TABFJW(ISTART + 121))*WDIF
     *                        + TABFJW(ISTART)
         D2WAL = WVAL + WVAL
         REXPW = DEXP(-WVAL)
         DENOM = FLOAT(JMAX0 + JMAX0 + 1)
         DO 100 J = JMAX0,1,-1
            DENOM = DENOM - TWO
            FJW(J - 1) = (D2WAL*FJW(J) + REXPW)/DENOM
  100    CONTINUE
      ELSE IF (IPOINT .LE. 20*JMAX0 + 360) THEN
         RWVAL = ONE/WVAL
         REXPW = DEXP(-WVAL)
         IRANGE = IPOINT/30 - 3
         GO TO (1215,1518,1824,1824,2430,2430), IRANGE
            FJW(0) = SQRPIH/SQRT(WVAL)
            GO TO 1000
 1215    CONTINUE
            GVAL = GFAC30
     *             + RWVAL*(GFAC31 + RWVAL*(GFAC32 + RWVAL*GFAC33))
            FJW(0) = SQRPIH/SQRT(WVAL) - REXPW*GVAL*RWVAL
            GO TO 1000
 1518    CONTINUE
            GVAL = GFAC20 + RWVAL*(GFAC21 + RWVAL*GFAC22)
            FJW(0) = SQRPIH/SQRT(WVAL) - REXPW*GVAL*RWVAL
            GO TO 1000
 1824    CONTINUE
            GVAL = GFAC10 + RWVAL*GFAC11
            FJW(0) = SQRPIH/SQRT(WVAL) - REXPW*GVAL*RWVAL
            GO TO 1000
 2430    CONTINUE
            GVAL = GFAC00
            FJW(0) = SQRPIH/SQRT(WVAL) - REXPW*GVAL*RWVAL
 1000    CONTINUE
         FACTOR = HALF*RWVAL
         TERM = FACTOR*REXPW
         DO 200 J = 1,JMAX0
            FJW(J) = FACTOR*FJW(J - 1) - TERM
            FACTOR = RWVAL + FACTOR
  200    CONTINUE
      ELSE
         FJW(0) = SQRPIH/SQRT(WVAL)
         RWVAL = ONE/WVAL
         FACTOR = HALF*RWVAL
         DO 300 J = 1,JMAX0
            FJW(J) = FACTOR*FJW(J-1)
            FACTOR = RWVAL + FACTOR
  300    CONTINUE
      END IF
      RETURN
C
C     ***** Tabulation of incomplete gamma function *****
C
      ENTRY GAMTAB(JMX)
C
C     For J = JMX a power series expansion is used, see for
C     example Eq.(39) given by V. Saunders in "Computational
C     Techniques in Quantum Chemistry and Molecular Physics",
C     Reidel 1975.  For J < JMX the values are calculated
C     using downward recursion in J.
C
C
      JMAX = JMX + 6
      IF (JMAX .GT. MAXJ) THEN
         WRITE (LUPRI,'(//A,I5,A,I3)')
     *      ' GAMTAB ERROR: JMAX =',JMAX,', which is greater than',MAXJ
         STOP
      END IF
      MAXJ0 = JMAX
C
C     WVAL = 0.0
C
      IADR = 1
      DENOM = ONE
      DO 700 J = 0,JMAX
         TABFJW(IADR) = ONE/DENOM
         IADR = IADR + 121
         DENOM = DENOM + TWO
  700 CONTINUE
C
C     WVAL = 0.1, 0.2, 0.3,... 12.0
C
      IADR = IADR - 121
      D2MAX1 = FLOAT(JMAX + JMAX + 1)
      R2MAX1 = ONE/D2MAX1
      DO 800 IPOINT = 1,120
         WVAL = TENTH*FLOAT(IPOINT)
         D2WAL = WVAL + WVAL
         IADR = IADR + 1
         TERM = R2MAX1
         SUM = TERM
         DENOM = D2MAX1
         DO 810 IORDER = 2,200
            DENOM = DENOM + TWO
            TERM = TERM*D2WAL/DENOM
            SUM = SUM + TERM
            IF (TERM .LE. 1.0D-15) GO TO 820
  810    CONTINUE
  820    CONTINUE
         REXPW = DEXP(-WVAL)
         TABFJW(IADR) = REXPW*SUM
         DENOM = D2MAX1
         JADR = IADR
         DO 830 J = 1,JMAX
            DENOM = DENOM - TWO
            TABFJW(JADR - 121) = (TABFJW(JADR)*D2WAL + REXPW)/DENOM
            JADR = JADR - 121
  830    CONTINUE
  800 CONTINUE
      RETURN
      END
