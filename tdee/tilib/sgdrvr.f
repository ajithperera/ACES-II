      SUBROUTINE SGDRVR (F, NEQN, Y, YP, T, TOUT, RELERR, ABSERR,
     $                   MAXNUM, IFLAG,
     $                   KODSTP, RMNSTP, RMXSTP, RSTPSM, RSTPSQ, ISTATS,
     $                   YY, WT, PHI, P, YPOUT)
C
C
C     Description of routine.
C
C     Differential equation solver routines for
C     Adams-Bashforth-Moulton predictor-ocorrector methods.
C
C     Subroutine SGDRVR integrates a system of up to NEQN first order
C     ordinary differential equations of the form
C
C     DY(I)/DT = F(T,Y(1),Y(2),...,Y(NEQN))
C
C     with Y(I) given at T.
C     The subroutine integrates from T to TOUT.  On return the
C     parameters in the call list are initialized for continuing the
C     integration.  The user has only to define a new value TOUT
C     and call SGDRVR again.
C
C     SGDRVR calls two routines, the integrator SGSTEP and the
C     interpolation routine SGNTRP. SGSTEP uses a modified divided
C     difference form of the Adams PECE
C     (Predict-Evaluate-Correct-Evaluate) formulas
C     and local extrapolation.  It adjusts the order
C     and step size to control the local error.  Normally each call
C     to SGSTEP advances the solution one step in the direction of TOUT.
C     for reasons of efficiency SGDRVR integrates beyond TOUT internally
C     though never beyond T + 10*(TOUT-T), and calls SGNTRP to
C     interpolate the solution at TOUT. An option is provided to stop
C     the integration at TOUT but it should be used only if it is
C     impossible to continue the integration beyond TOUT .
C
C     This code is completely explained and documented in the text,
C     Computer solution of ordinary differential equations:  the initial
C     value problem by L. F. Shampine and M. K. Gordon (Freeman, 1975)
C
C     The parameters for SGDRVR are:
C
C     F -- subroutine F(T,Y,YP) to evaluate derivatives YP(I)=dY(I)/dT
C
C     NEQN -- number of equations to be integrated
C
C     Y(*) -- solution vector at T
C
C     T -- independent variable
C
C     TOUT -- point at which solution is desired
C
C     RELERR,ABSERR -- relative and absolute error tolerances for local
C            error test.  at each step the code requires
C     
C            ABS(local error) .LE. ABS(Y)*RELERR + ABSERR
C     
C            for each component of the local error and solution vectors
C
C     MAXNUM -- maximum number of steps allowed in one call to SGDRVR
C
C     IFLAG -- indicaes status of integration
C
C     First call to SGDRVR
C
C     The user must provide cstorage in his calling program for the
C     array in the call list Y(NEQN), declare F in an external
C     statement, supply the subroutine F(T,Y,YP) to evaluate
C
C     dY(I)/dT = YP(I) = F(T,Y(1),Y(2),...,Y(NEQN))
C
C     and initialize the parameters:
C
C     NEQN -- number of equations to te integrated
C
C     Y(*) -- vector of initial conditions
C
C     T -- starting point of integration
C
C     TOUT -- point at which solution is desired
C
C     RELERR,ABSERR -- relative and absolute local error tolerances
C
C     IFLAG -- +1,-1.  indicator to initialize the code.  Normal input
C          is +1.  the user should set IFLAG=-1 only if it is
C          impossible to continue the integration beyond TOUT .
C
C     All parameters except F , NEQN and TOUT may be altered by the
C     code on output so must be variables in the calling program.
C
C     Output from SGDRVR
C
C     NEQN -- unchanged
C
C     Y(*) -- solution at T
C
C     T -- last point reached in integration.  Normal return has
C          T = TOUT
C
C     TOUT -- unchaged
C
C     RELERR,ABSERR -- normal return has tolerances unchanged.  IFLAG=3
C           signals tolerances increased
C
C     IFLAG = 2 -- normal return.  integration reached TOUT
C
C           = 3 -- integration did not reach tout because error
C                  tolerances too small.  RELERR, ABSERR increased
C                  appropriately for continuing
C
C           = 4 -- integration did not reach TOUT because more than
C                  maximum steps needed
C
C           = 5 -- integration did not reach TOUT because equations
C                  appear to be stiff
C
C           = 6 -- invalid input parameters (fatal error)
C                  the value of IFLAG is returned negative when the
C                  input value is negative and the integration does
C                  not reach TOUT, i.e., -3, -4, -5.
C
C     Subsequent calls to SGDRVR
C
C     Subroutine SGDRVR returns with all information needed to continue
C     the integration.  If the integration reached TOUT , the user need
C     only define a new TOUT and call again.  If the integration did not
C     reach TOUT and the user wants to continue, he just calls again.
C     The output value of IFLAG is the appropriate input value for
C     subsequent calls.  the only situation in which it should be
C     altered is to stop the integration internally at the new TOUT,
C     i.e., change output IFLAG=2 to input IFLAG=-2 .  Error tolerances
C     may be changed by the user before continuing.  All other
C     parameters must remain unchanged.
C
C     Original code: copyright Shampine and Gordon 1975.
C     Modifications: copyright Erik Deumens and QTP, 1990, 1996.
C
      IMPLICIT NONE
C
C     The only machine dependent constant is based on the machine unit
C     roundoff error  U which is the smallest positive number such that
C     1.0+U .GT. 1.0 .  U must be calculated and FOURU=4.0*U  inserted
C     in the following data statement before using DE.  The routine
C     MACHIN calculates U as follow:
C
C     .     HALFU = HALF
C     .  50 TEMP1 = ONE + HALFU
C     .     IF(TEMP1.LE.ONE) GO TO 100
C     .     HALFU = HALF*HALFU
C     .     GO TO 50
C     . 100 U     = TWO*HALFU
C
C     FOURU and TWOU=2.0*U must also be inserted in subroutine SGSTEP
C     before calling SGDRVR.
C     These machine dependent values are chosen for the Cray.
C     Sun and RS/6000 could have smaller numbers but the
C     difference is not important.
C
      DOUBLE PRECISION FOURU, ZERO, ONE, TEN  
      PARAMETER (FOURU = 2.8D-14, ZERO = 0.D0, ONE = 1.D0,
     $     TEN = 10.D0)
C
C     Arguments
C
      EXTERNAL F
      INTEGER NEQN, ISTATS, MAXNUM, IFLAG
      DOUBLE PRECISION
     $     Y(NEQN), YP(NEQN), T, TOUT, RELERR, ABSERR
      DOUBLE PRECISION
     $     YY(NEQN), WT(NEQN), PHI(NEQN,16), P(NEQN), YPOUT(NEQN)
      INTEGER KODSTP
      DOUBLE PRECISION RMNSTP, RMXSTP, RSTPSM, RSTPSQ
C
C     Local variables
C
      DOUBLE PRECISION TOLD, H, HOLD, X, TEND, PSI(12)
      DOUBLE PRECISION ABSEPS, RELEPS, EPS, DEL, ABSDEL   
      LOGICAL START,CRASH,STIFF
      INTEGER KOLD, L, NOSTEP, ISN, ISNOLD, KLE4, K
      SAVE TOLD, START, CRASH, STIFF, PSI, ISNOLD, H, X, K
C-----------------------------------------------------------------------
C
C     Test for improper parameters
C
      IF ( NEQN.LT.1 .OR. T.EQ.TOUT
     $     .OR. RELERR.LT.ZERO .OR. ABSERR.LT.ZERO ) THEN
        IFLAG = 6
        RETURN
      ELSE
        EPS = MAX(RELERR,ABSERR)
      END IF
C
      IF ( IFLAG .EQ. 0 ) THEN
        IFLAG = 6
        RETURN
      ELSE
        ISN = ISIGN(1,IFLAG)
        IFLAG = IABS(IFLAG)
      END IF
C
      IF ( (T.NE.TOLD .AND. IFLAG.NE.1)
     $              .OR.
     $     (IFLAG.LT.1 .OR. IFLAG.GT.5) ) THEN
        IFLAG = 6
        RETURN
      END IF
C
C     On each call set interval of integration and counter for number
C     of steps.  Adjust input error tolerances to define weight vector
C     for subroutine  SGSTEP.
C
      DEL = TOUT - T
      ABSDEL = ABS(DEL)
      TEND = T + TEN*DEL
      IF(ISN .LT. 0) TEND = TOUT
      NOSTEP = 0
      KLE4 = 0
      STIFF = .FALSE.
      RELEPS = RELERR/EPS
      ABSEPS = ABSERR/EPS
C
      IF ( IFLAG.EQ.1 .OR. ISNOLD.LT.0 .OR. ABSDEL.EQ.ZERO) THEN
C
C       On start and restart also set work variables X and YY(*),
C       store the direction of integration and initialize the step size
C
        START = .TRUE.
        X = T
        CALL DCOPY(NEQN,Y,1,YY,1)
        H = SIGN(MAX(ABS(TOUT-X),FOURU*ABS(X)),TOUT-X)
      END IF
C
C     While ( .NOT. CRASH ) do ...
C
   50 CONTINUE
C
      IF ( .NOT. (ABS(X-T) .LT. ABSDEL) ) THEN
C
C       If already past output point, interpolate and return
C
        CALL SGNTRP (X,YY,TOUT,Y,YPOUT,NEQN,KOLD,PHI,PSI)
        IFLAG = 2
        T = TOUT
        TOLD = T
        ISNOLD = ISN
        RETURN
      END IF
C
C     If cannot go past output point and sufficiently close,
C     extraplolate and return
C
      IF ( .NOT. (ISN.GT.0  .OR.  ABS(TOUT-X).GE.FOURU*ABS(X)) ) THEN
        H = TOUT - X
        CALL F(X,YY,YP,ISTATS)
        IF (MOD(ISTATS,2) .NE. 1) THEN
          RETURN
        ELSE
          CALL DCOPY(NEQN,YY,1,Y,1)
          CALL DAXPY(NEQN,H,YP,1,Y,1)
          IFLAG = 2
          T = TOUT
          TOLD = T
          ISNOLD = ISN
          RETURN
        END IF
      END IF
C
C     Test for too much work
C
      IF ( NOSTEP .GE. MAXNUM) THEN
        IFLAG = ISN*4
        IF ( STIFF ) IFLAG = ISN*5
        CALL DCOPY(NEQN,YY,1,Y,1)
        T = X
        TOLD = T
        ISNOLD = 1
        RETURN
      END IF
C
C     Limit step size, set weight vector and take a step
C
      H = SIGN(MIN(ABS(H),ABS(TEND-X)),H)
      DO 110 L = 1,NEQN
        WT(L) = RELEPS*ABS(YY(L)) + ABSEPS
  110 CONTINUE
      CALL SGSTEP ( X,YY,F,ISTATS,NEQN,H,EPS,WT,START,
     $              HOLD,K,KOLD,CRASH,PHI,P,YP,PSI    )
      IF (MOD(ISTATS,2) .NE. 1) RETURN
C
C     Test for tolerances too small
C
      IF ( CRASH ) THEN
        IFLAG = ISN*3
        RELERR = EPS*RELEPS
        ABSERR = EPS*ABSEPS
        CALL DCOPY(NEQN,YY,1,Y,1)
        T = X
        TOLD = T
        ISNOLD = 1
        RETURN
      ELSE
C
C       Augment counter on work and test for stiffness
C
        NOSTEP = NOSTEP + 1
        KODSTP = KODSTP + 1
        RMNSTP = MIN(RMNSTP,HOLD)
        RMXSTP = MAX(RMXSTP,HOLD)
        RSTPSM = RSTPSM + HOLD
        RSTPSQ = RSTPSQ + HOLD**2
        KLE4 = KLE4 + 1
        IF(KOLD .GT. 4) KLE4 = 0
        IF(KLE4 .GE. 50) STIFF = .TRUE.
        GO TO 50
      END IF
      END
