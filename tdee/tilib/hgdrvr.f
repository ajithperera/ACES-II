      SUBROUTINE HGDRVR (DIFFUN, JAC, N, T0, H0, Y0, TOUT, EPS,
     $     MF, INDEX, ML, MU,
     $     KODSTP, RMNSTP, RMXSTP, RSTPSM, RSTPSQ, ISTATS,
     $     Y, YMAX, ERROR, SAVE1, SAVE2, PW, IPIV)
C
C     Description of routine.
C
C     This is the March 19, 1975 version of
C     GEAR-B, a package for the solution of the initial value
C     problem for systems of ordinary differential equations,
C     
C     dy/dt = f(y,t),    y = (y(1),y(2),...,y(n))
C     .
C     GEAR-B is a variant of the GEAR package to be used when
C     the jacobian matrix df/dy has banded or nearly banded form.
C     subroutine HGDRVR is a driver routine for the GEAR-B package.
C     GEAR-B uses Adams-Bashforth-Moulton predictor-corrector methods
C     for non-stiff problems and Gear's backward differentiation
C     formulas (BDF) for stiff problems.
C
C     References:
C
C     1.  A. C. Hindmarsh, Gear: Grdinary differential equation
C             system solver, ucid-30001 rev. 3, Lawrence Livermore
C             Laboratory, p.o.box 808, Livermore, CA 94550, Dec. 1974.
C
C     2.  A. C. Hindmarsh, Gearb: Solution of ordinary
C             differential equations having banded jacobian,
C             ucid-30059 rev. 1, l.l.l., March 1975.
C
C     HGDRVR is to be called once for each output value of t, and
C     in turn makes repeated calls to the core integrator, HGSTEP.
C
C     The input parameters are:
C
C     N     =  The number of first-order differential equations.
C               n can be reduced, but never increased, during a problem.
C
C     T0    =  The initial value of t, the independent variable
C               (used only on first call).
C
C     H0    =  The next step size in t (used for input only on the
C               first call).
C
C     Y0    =  A vector of length n containing the initial values of
C               y (used for input only on first call).
C
C     TOUT  =  The value of t at which output is desired next.
C               integration will normally go slightly beyond tout
C               and the package will interpolate to t = tout.
C
C     EPS   =  The relative error bound  (used only on the
C               first call, unless index = -1).  Single step error
C               estimates divided by YMAX(I) will be kept less than
C               EPS in maximum. The vector YMAX of weights
C               is computed in HGDRVR.  YMAX(I) is ABS(Y(I))+ONE.
C               This defines the absolute error for small components
C               to be EPS.
C
C     MF    =  The method flag  (used only on first call, unless
C               INDEX = -1).  Allowed values are 10, 11, 12, 13,
C               20, 21, 22, 23.  MF has two decimal digits, METH
C               and MITER  (MF = 10*METH + MITER).
C               METH is the basic method indicator:
C
C                 METH = 1  means the Adams methods.
C
C                 METH = 2  means the backward differentiation
C                           formulas (BDF), or stiff methods of Gear.
C
C               MITER is the iteration method indicator:
C
C                 MITER = 0 means functional iteration (no partial
C                           derivatives needed).
C
C                 MITER = 1 means chord method with analytic jacobian.
C                           for this user supplies subroutine
C                           JAC  (see description below).
C
C                 MITER = 2 means chord method with jacobian calculated
C                           internally by finite differences.
C
C                 MITER = 3 means chord method with jacobian replaced
C                           by a diagonal approximation based on a
C                           directional derivative.
C
C     INDEX =  Integer used on input to indicate type of call,
C              with the following values and meanings:
C
C                  1    this is the first call for this problem.
C
C                  0    this is not the first call for this problem,
C                       and integration is to continue.
C
C                 -1    this is not the first call for the problem,
C                       and the user has reset N, EPS, and/or MF.
C
C                  2    same as 0 except that TOUT is to be hit
C                       exactly (no interpolation is done).
C                       assumes TOUT .GE. the current t.
C
C                  3    same as 0 except control returns to calling
C                       program after one step.  TOUT is ignored.
C
C               Since the normal output value of INDEX is 0,
C               it need not be reset for normal continuation.
C
C     ML,MU =  The widths of the lower and upper parts, respectively,
C               of the band in the jacobian matrix, not counting the
C               main diagonal.  The full bandwidth is ML + MU + 1.
C
C     After the initial call, if a normal return occurred and a normal
C     continuation is desired, simply reset TOUT and call again.
C     All other parameters will be ready for the next call.
C     A change of parameters with INDEX = -1 can be made after
C     either a successful or an unsuccessful return.
C
C     the output parameters are:
C
C     H0    =  The step size H used last, whether successfully or not.
C
C     Y0    =  The computed values of y at t = tout.
C
C     TOUT  =  The output value of t.  If integration was successful,
C               and the input value of index was not 3, TOUT is
C               unchanged from its input value.  Otherwise, tout
C               is the current value of t to which integration
C               has been completed.
C
C     INDEX =  Integer used on output to indicate results,
C               with the following values and meanings:
C
C          0    Integration was completed to TOUT or beyond.
C
C         -1    The integration was halted after failing to pass the
C               error test even after reducing h by a factor of
C               1.E10 from its initial value.
C
C         -2    After some initial success, the integration was
C               halted either by repeated error test failures or by
C               a test on EPS.  Too much accuracy has been requested.
C
C         -3    The integration was halted after failing to achieve
C               corrector convergence even after reducing h by a
C               factor of 1.E10 from its initial value.
C
C         -4    Immediate halt because of illegal values of input
C               parameters.  See printed message.
C
C         -5    Index was -1 on input, but the desired changes of
C               parameters were not implemented because tout
C               was not beyond t.  Interpolation to t = TOUT was
C               performed as on a normal return.  to try again,
C               simply call again with index = -1 and a new tout.
C
C     In addition to HGDRVR, the following routines are provided in
C     the package:
C
C     HGNTRP(TOUT,Y,N0,Y0)  Interpolates to get the output values
C                           at t = TOUT, from the data in the y array.
C
C     HGSTEP(Y,N0)  is the core integrator routine.  It performs a
C                           single step and associated error control.
C
C     DIFFUN and JAC passed as argument and DIFFUN changed to have
C     uniform argument list for integrators
C
C     call is now : HGSTEP(DIFFUN,JAC,ISTATS,Y,N0,
C                          YMAX, ERROR, SAVE1, SAVE2, PW, IPIV)
C     (ED QTP 30-Dec-1986)
C
C     HGCOEF(METH,NQ,EL,TQ,MAXDER)  sets coefficients for use in
C                           the core integrator.
C
C     HGPSET(Y,N0,CON,MITER,IER)  computes and processes the jacobian
C                           matrix J = df/dy.
C
C     call is now : HGPSET(DIFFUN,JAC,ISTATS,Y,N0,CON,MITER,IER,
C                         YMAX, ERROR, SAVE1, SAVE2, PW, IPIV)
C     (ED QTP 30-Dec-1986)
C
C     HGDEC(N0,N,ML,MU,B,IP,IER)  and  HGSOL(N0,N,ML,MU,B,C,X,IP)  are
C                           used to solve banded linear systems.
C
C     Note: HGPSET, HGDEC, and HGSOL are called only if MITER = 1 or 2.
C
C     The following routines are to be supplied by the user:
C
C     DIFFUN (N,T,Y,YDOT)   computes the function YDOT = f(y,t), the
C                           right-hand side of the o.d.e.
C                           here y and ydot are vectors of length N.
C
C     call is now : DIFFUN (T,Y,YDOT,ISTATS)
C     (ED QTP 30-Dec-1986)
C
C     JAC(N,T,Y,PD,N0,ML,MU)  computes the N by N jacobian matrix of
C                           partial derivatives, and stores it in PD
C                           as an N0 by ML+MU+1 array.  PD(I,J-I+ML+1)
C                           is to be set to the partial derivative of
C                           YDOT(I) with respect to Y(J).  Called only
C                           if MITER = 1.  Otherwise a dummy routine
C                           can be substituted.
C
C     Arrays are passed as arguments to make them dynamic
C     (ED QTP 31-Dec-1986)
C     
C     Dimension Y(N,13)   if METH <> 2
C               Y(N, 6)   if METH = 2
C
C     Dimension YMAX(N), ERROR(N), SAVE1(N), SAVE2(N),
C               PW(N*(2*ML+MU+1))   if MITER = 1 or 2
C               PW(N)               if MITER = 3
C               PW(1)               if MITER = 0
C               IPIV(N)   if MITER = 1 or 2
C               IPIV(1)   if MITER = 0 or 3
C     
C
C     The dimensions in the following declarations are set for a
C     maximum of 100 equations and for 2*ML+MU .le. 30.  If these limits
C     are to be exceeded, the dimensions should be increased accordingly
C     the dimension of PW below must be at least N*(2*ML+MU+1) if MITER
C     1 or 2, but can be reduced to N if MITER = 3, or to 1 if MITER = 0
C     The dimensions of YMAX, ERROR, SAVE1, SAVE2, IPIV, and the first
C     dimension of Y should all be at least N.  The column length of
C     the Y array as used elsewhere is N0, not 100.  The row length of Y
C     can be reduced from 13 to 6 if METH = 2.
C     The IPIV array is used only if MITER is 1 or 2.
C
C     The common block GEAR9 can be accessed externally by the user
C     if desired.  It contains the step size last used (successfully),
C     the order last used (successfully), the number of steps taken
C     so far, the number of f evaluations (DIFFUN calls) so far,
C     and the number of jacobian evaluations so far.
C
C     Original code: copyright A. C. Hindmarsh and LLL, 1975.
C     Modifications: copyright Erik Deumens and QTP, 1986, 1996.
C
      IMPLICIT NONE

      INTEGER LOUT
      PARAMETER (LOUT = 6)
      DOUBLE PRECISION ONE
      PARAMETER (ONE = 1.D0)
C
C     Arguments
C
      EXTERNAL DIFFUN, JAC
      INTEGER N, MF, INDEX, ML, MU, ISTATS
      DOUBLE PRECISION Y0(N), Y(N,13)
      DOUBLE PRECISION T0, H0, TOUT, EPS
      DOUBLE PRECISION
     $     YMAX(N), ERROR(N), SAVE1(N), SAVE2(N), PW(N*(2*ML+MU+1))
      INTEGER IPIV(N)
      INTEGER KODSTP
      DOUBLE PRECISION RMNSTP, RMXSTP, RSTPSM, RSTPSQ
C
C     Global variables for HG routines.
C
      DOUBLE PRECISION T, H, HMIN, HMAX, EPSC, UROUND
      INTEGER NC, MFC, KFLAG, JSTART
      COMMON /GEAR1/ T,H,HMIN,HMAX,EPSC,UROUND,NC,MFC,KFLAG,JSTART
      DOUBLE PRECISION EPSJ
      INTEGER MLC, MUC, MW, NM1, N0ML, N0W
      COMMON /GEAR8/ EPSJ,MLC,MUC,MW,NM1,N0ML,N0W
      DOUBLE PRECISION HUSED
      INTEGER NQUSED, NSTEP, NFE, NJE
      COMMON /GEAR9/ HUSED,NQUSED,NSTEP,NFE,NJE
C
C     Local variables
C
      DOUBLE PRECISION TOUTP, AYI
C     Erik Apr 1996
C*CALL REAL
C     $     D
      INTEGER I, KGO, NHCUT, N0
      SAVE N0
C-----------------------------------------------------------------------
C     UROUND =  The unit roundoff of the machine, i.e. the smallest
C              positive U such that 1. + U .NE. 1. on the machine.
C
      UROUND = 1.D-14
      IF (INDEX .EQ. 0) GO TO 20
      IF (INDEX .EQ. 2) GO TO 25
      IF (INDEX .EQ. -1) GO TO 30
      IF (INDEX .EQ. 3) GO TO 40
      IF (INDEX .NE. 1) GO TO 430
      IF (EPS .LE. UROUND) GO TO 400
      IF (N .LE. 0) GO TO 410
      IF ((T0-TOUT)*H0 .GE. 0.D0) GO TO 420
C
C     If initial values of YMAX other than those set below are desired,
C     they should be set here.  All YMAX(I) must be positive.
C     If values for HMIN or HMAX, the bounds on ABS(H), other than
C     those below are desired, they should be set below.
C
      DO 10 I = 1,N
C       ED Apr 1996 >>>
C       Make YMAX a scale for local error.
C       YMAX(I) = MAX(ABS(Y0(I)),ONE)
        YMAX(I) = ABS(Y0(I)) + ONE
C       ED Apr 1996 <<<
 10     Y(I,1) = Y0(I)
      NC = N
      T = T0
      H = H0
      IF ((T+H) .EQ. T) WRITE(LOUT,15)
 15   FORMAT(1X,' @HGDRVR-W,',' Warning: T + H = T on next step.')
      HMIN = ABS(H0)
      HMAX = ABS(T0-TOUT)*10.D0
      EPSC = EPS
      MFC = MF
      JSTART = 0
      N0 = N
      EPSJ = SQRT(UROUND)
      MLC = ML
      MUC = MU
      MW = ML + MU + 1
      NM1 = N0 - 1
      N0ML = N0*ML
      N0W = N0*MW
      NHCUT = 0
      GO TO 50
C
C     TOUTP is the previous value of TOUT for use in HMAX.
 20   HMAX = ABS(TOUT-TOUTP)
      GO TO 80
C
 25   HMAX = ABS(TOUT-TOUTP)
      IF ((T-TOUT)*H .GE. 0.D0) GO TO 500
      GO TO 85
C
 30   IF ((T-TOUT)*H .GE. 0.D0) GO TO 440
      JSTART = -1
      NC = N
      EPSC = EPS
      MFC = MF
C
 40   IF ((T+H) .EQ. T) WRITE(LOUT,15)
C
 50   CALL HGSTEP (DIFFUN, JAC, ISTATS, Y, N0,
     $      YMAX, ERROR, SAVE1, SAVE2, PW, IPIV)
      IF (MOD(ISTATS,2) .NE. 1) RETURN
      KODSTP = KODSTP + 1
      RMNSTP = MIN(RMNSTP,HUSED)
      RMXSTP = MAX(RMXSTP,HUSED)
      RSTPSM = RSTPSM + HUSED
      RSTPSQ = RSTPSQ + HUSED**2
C
      KGO = 1 - KFLAG
      GO TO (60, 100, 200, 300), KGO
C     KFLAG  =   0,  -1,  -2,  -3
C
 60   CONTINUE
C
C     Normal return from integrator.
C
C     The weights YMAX(I) are updated.  If different values are desired,
C     they should be set here.  A test is made for EPS being too small
C     for the machine precision.
C
C     Any other tests or calculations that are required after every
C     step should be inserted here.
C
C     If INDEX = 3, Y0 is set to the current Y values on return.
C     If INDEX = 2, H is controlled to hit TOUT (within roundoff
C     error), and then the current Y values are put in Y0 on return.
C     For any other value of INDEX, control returns to the integrator
C     unless TOUT has been reached.  Then interpolated values of Y are
C     computed and stored in Y0 on return.
C     If interpolation is not desired, the call to HGNTRP should be
C     removed and control transferred to statement 500 instead of 520.
C
C     Test whether requested accuaracy is attainable.
C
C     D = 0.D0
      DO 70 I = 1,N
        AYI = ABS(Y(I,1))
C       ED Apr 96 >>>
C       YMAX(I) = MAX(YMAX(I), AYI)
C       D = D + (AYI/YMAX(I))**2
        YMAX(I) = AYI + ONE
   70 CONTINUE
C     D = D*(UROUND/EPS)**2
C     IF (D .GT. FLOAT(N)) GO TO 250
      IF (UROUND/EPS .GT. 1.D0) GO TO 250
C     ED Apr 96 <<<
      IF (INDEX .EQ. 3) GO TO 500
      IF (INDEX .EQ. 2) GO TO 85
 80   IF ((T-TOUT)*H .LT. 0.D0) GO TO 40
      CALL HGNTRP (TOUT, Y, N0, Y0)
      GO TO 520
 85   IF (((T+H)-TOUT)*H .LE. 0.D0) GO TO 40
      IF (ABS(T-TOUT) .LE. 100.D0*UROUND*HMAX) GO TO 500
      IF ((T-TOUT)*H .GE. 0.D0) GO TO 500
      H = (TOUT - T)*(1.D0 - 4.D0*UROUND)
      JSTART = -1
      GO TO 40
C
C     On an error return from integrator, an immediate return occurs if
C     KFLAG = -2, and recovery attempts are made otherwise.
C     To recover, H and HMIN are reduced by a factor of .1 up to 10
C     times before giving up.
C
 100  WRITE (LOUT,105) T
 105  FORMAT(1X,' @HGDRVR-W',' KFLAG = -1 from integrator at T = ',
     $     E16.8/10X,' Error test failed with ABS(H) = HMIN')
 110  IF (NHCUT .EQ. 10) GO TO 150
      NHCUT = NHCUT + 1
      HMIN = .1D0*HMIN
      H = .1D0*H
      WRITE (LOUT,115) H
 115  FORMAT(1X,' @HGDRVR-I',' H has been reduced to ',E16.8,
     $     ' and step will be retried.')
      JSTART = -1
      GO TO 40
C
 150  WRITE (LOUT,155)
 155  FORMAT(1X,' @HGDRVR-F',
     $     ' Problem appears unsolvable with given input.')
      GO TO 500
C
 200  WRITE (LOUT,205) T,H
 205  FORMAT(1X,' @HGDRVR-F',
     $     ' KFLAG = -2 from integrator at T = ',E16.8,' H = ',E16.8/
     $     10X,' The requested error is smaller than can be handled.')
      GO TO 500
C
 250  WRITE (LOUT,255) T
 255  FORMAT(1X,' @HGDRVR-F',
     $     ' Integration halted by driver at T = ',E16.8/10X,
     $     ' EPS too small to be attained for the machine precision.')
      KFLAG = -2
      GO TO 500
C
 300  WRITE (LOUT,305) T
  305 FORMAT(1X,' @HGDRVR-W',' KFLAG = -3 from integrator at T = ',
     $     E16.8/10X,
     $     ' Corrector convergence could not be achieved.')
      GO TO 110
C
 400  WRITE (LOUT,405) UROUND
  405 FORMAT(1X,' @HGDRVR-F',' Illegal input:  EPS .LE. ',E16.8)
      INDEX = -4
      RETURN
C
 410  WRITE (LOUT,415)
 415  FORMAT(1X,' @HGDRVR-F',' Illegal input:  N .LE. 0')
      INDEX = -4
      RETURN
C
 420  WRITE (LOUT,425)
 425  FORMAT(1X,' @HGDRVR-F',' Illegal input:  (T0-TOUT)*H .GE. 0.')
      INDEX = -4
      RETURN
C
 430  WRITE (LOUT,435) INDEX
 435  FORMAT(1X,' @HGDRVR-F',' Illegal input:  INDEX = ',I5)
      INDEX = -4
      RETURN
C
 440  WRITE(LOUT,445) T,TOUT,H
 445  FORMAT(1X,' @HGDRVR-I',
     $     ' INDEX = -1 on input with (T-TOUT)*H .GE. 0.'/
     $     10X,' T = ',E16.8,' TOUT = ',E16.8,' H = ',E16.8/
     $     10X,' Interpolation was done as on normal return.'/
     $     10X,' Desired parameter changes were not made.')
      CALL HGNTRP (TOUT, Y, N0, Y0)
      INDEX = -5
      RETURN
C
 500  TOUT = T
      DO 510 I = 1,N
 510    Y0(I) = Y(I,1)
 520  INDEX = KFLAG
      TOUTP = TOUT
      H0 = HUSED
      IF (KFLAG .NE. 0) H0 = H
      RETURN
      END
