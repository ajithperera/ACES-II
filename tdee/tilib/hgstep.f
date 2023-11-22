      SUBROUTINE HGSTEP (DIFFUN, JAC, ISTATS, Y, N0,
     $      YMAX, ERROR, SAVE1, SAVE2, PW, IPIV)
C
C     Description of routine.
C     HGSTEP performs one step of the integration of an initial value
C     problem for a system of ordinary differential equations.
C     HGSTEP is a version for banded form of the jacobian matrix.
C     communication with HGSTEP is done with the following variables.
C
C     Y       An N0 by LMAX array containing the dependent variables
C             and their scaled derivatives.  LMAX is 13 for the Adams
C             methods and 6 for the Gear methods.  LMAX - 1 = MAXDER
C             is the maximum order available.  See subroutine HGCOEF.
C             Y(I,J+1) contains the J-th derivative of Y(I), scaled by
C             H**J/FACTORIAL(J)  (J = 0,1,...,NQ).
C
C     N0      A constant integer .GE. N, used for dimensioning purposes.
C
C     T       The independent variable. t is updated on each step taken.
C
C     H       The step size to be attempted on the next step.
C             H is altered by the error control algorithm during the
C             problem.  H can be either positive or negative, but its
C             sign must remain constant throughout the problem.
C
C     HMIN, HMAX   The minimum and maximum absolute value of the step
C             size to be used for the step.  These may be changed at
C             any time, but will not take effect until the next H
C             change.
C
C     EPS     The relative error bound.  See description in DRIVEB.
C
C     UROUND  The unit roundoff of the machine.
C
C     N       The number of first-order differential equations.
C
C     MF      The method flag.  See description in DRIVEB.
C
C     KFLAG   A completion code with the following meanings.
C
C                      0  the step was succesful.
C
C                     -1  the requested error could not be achieved
C                           with ABS(H) = HMIN.
C
C                     -2  the requested error is smaller than can
C                           be handled for this problem.
C
C                     -3  corrector convergence could not be
C                           achieved for ABS(H) = HMIN.
C
C             On a return with KFLAG negative, the values of T and
C             the Y array are as of the beginning of the last
C             step, and H is the last step size attempted.
C
C     JSTART  An integer used on input and output.
C             On input, it has the following values and meanings.
C
C                      0  perform the first step.
C
C                  .GT.0  take a new step continuing from the last.
C
C                  .LT.0  take the next step with a new value of
C                           H, EPS, N, and/or MF.
C
C             On exit, JSTART is NQ, the current order of the method.
C
C     YMAX    An array of N elements with which the estimated local
C             errors in Y are compared.
C
C     ERROR   An array of N elements.  ERROR(I)/TQ(2) is the estimated
C             one-step error in Y(I).
C
C     SAVE1, SAVE2  Two arrays of working storage,
C             each of length N.
C
C     PW      A block of locations used for partial derivatives if
C             MITER is not 0.  See description in DRIVEB.
C
C     IPIV    An integer array of length N used for pivot
C             information if MITER = 1 or 2.
C
C     ML,MU   The lower and upper half bandwidths, respectively, of
C             the jacobian.  See description in DRIVEB.
C
C     Original code: copyright A. C. Hindmarsh and LLL, 1975.
C     Modifications: copyright Erik Deumens and QTP, 1986, 1992.
C
C
      IMPLICIT NONE

C TVG discovered save dependence
      save

C
C     Global variables
C
      DOUBLE PRECISION T,H,HMIN,HMAX,EPS,UROUND
      INTEGER N,MF,KFLAG,JSTART
      COMMON /GEAR1/ T,H,HMIN,HMAX,EPS,UROUND,N,MF,KFLAG,JSTART
      DOUBLE PRECISION EPSJ
      INTEGER ML, MU, MW
      COMMON /GEAR8/ EPSJ,ML,MU,MW,NM1
      INTEGER NQUSED, NSTEP, NFE, NJE
      DOUBLE PRECISION HUSED
      COMMON /GEAR9/ HUSED,NQUSED,NSTEP,NFE,NJE
C
C     Arguments
C
      EXTERNAL DIFFUN, JAC
      INTEGER ISTATS
      INTEGER N0
      DOUBLE PRECISION Y(N0,13)
      DOUBLE PRECISION YMAX(N0), ERROR(N0), SAVE1(N0),
     $                 SAVE2(N0), PW(N0*(2*Ml+MU+1))
      INTEGER IPIV(N0)
C
C     Local variables
C
      DOUBLE PRECISION
     $   EL(13),TQ(4), TOLD, PR1, PR2, PR3, ENQ1, ENQ2, ENQ3, HL0, PHL0
      DOUBLE PRECISION
     $   D, D1, R, R0, R1, RH, RMAX, RC, CON, BND, E, EUP,
     $   EDN, CRATE, ZP9, EPSOLD, OLDL0, HOLD
C     ED Apr 1996
C*CALL REAL
C     $     FN
      INTEGER I, J, J1, J2, M, L, IDOUB, MFOLD, MEO, NM1, NEWQ
      INTEGER LMAX, IREDO, NQ
      SAVE NQ
      INTEGER NSTEPJ, NOLD, IER, METH, MAXDER, IWEVAL, MITER, MIO, IRET

C
      DATA EL(2)/1.D0/, OLDL0/1.D0/, ZP9 /.9D0/
C-----------------------------------------------------------------------
      KFLAG = 0
      TOLD = T
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .NE. 0) GO TO 120
C
C     On the first call, the order is set to 1 and the initial YDOT is
C     calculated.  RMAX is the maximum ratio by which H can be increased
C     in a single step.  It is initially 10,000. to compensate for the
C     small initial H, but then is normally equal to 10.  If a failure
C     occurs (in corrector convergence or error test), RMAX is set at 2
C     for the next increase.
C
      CALL DIFFUN (T, Y, SAVE1, ISTATS)
      IF (MOD(ISTATS,2) .NE. 1) RETURN
      DO 110 I = 1,N
 110    Y(I,2) = H*SAVE1(I)
      METH = MF/10
      MITER = MF - 10*METH
      NQ = 1
      L = 2
      IDOUB = 3
C     ED Apr 1996
C     RMAX = 1.D4
      RMAX = 10.D0
      RC = 0.D0
      CRATE = 1.D0
      EPSOLD = EPS
      HOLD = H
      MFOLD = MF
      NOLD = N
      NSTEP = 0
      NSTEPJ = 0
      NFE = 1
      NJE = 0
      IRET = 3
      GO TO 130
C
C     If the caller has changed METH, HGCOEF is called to set
C     the coefficients of the method.  If the caller has changed
C     N, EPS, or METH, the constants E, EDN, EUP, and BND must be reset.
C     E is a comparison for errors of the current order NQ. EUP is
C     to test for increasing the order, EDN for decreasing the order.
C     BND is used to test for convergence of the corrector iterates.
C     If the caller has changed H, Y must be rescaled.
C     Further changes in H for that many steps.
C     If H or METH has been changed, IDOUB is reset to L + 1 to prevent
C
 120  IF (MF .EQ. MFOLD) GO TO 150
      MEO = METH
      MIO = MITER
      METH = MF/10
      MITER = MF - 10*METH
      MFOLD = MF
      IF (MITER .NE. MIO) IWEVAL = MITER
      IF (METH .EQ. MEO) GO TO 150
      IDOUB = L + 1
      IRET = 1
 130  CALL HGCOEF (METH, NQ, EL, TQ, MAXDER)
      LMAX = MAXDER + 1
      RC = RC*EL(1)/OLDL0
      OLDL0 = EL(1)
 140  CONTINUE
C
C     ED Apr 1996 >>>
C     Use MAX norm instead of R.M.S. norm.
C      FN = FLOAT(N)
C      EDN = FN*(TQ(1)*EPS)**2
C      E   = FN*(TQ(2)*EPS)**2
C      EUP = FN*(TQ(3)*EPS)**2
C      BND = FN*(TQ(4)*EPS)**2
      EDN = TQ(1)*EPS
      E   = TQ(2)*EPS
      EUP = TQ(3)*EPS
      BND = TQ(4)*EPS
C     ED Apr 1996 <<<
      GO TO (160, 170, 200), IRET
 150  IF ((EPS .EQ. EPSOLD) .AND. (N .EQ. NOLD)) GO TO 160
      EPSOLD = EPS
      NOLD = N
      IRET = 1
      GO TO 140
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,HMAX/ABS(H),RMAX)
      R1 = 1.D0
      DO 180 J = 2,L
        R1 = R1*RH
        DO 180 I = 1,N
 180      Y(I,J) = Y(I,J)*R1
      H = H*RH
      RC = RC*RH
      IDOUB = L + 1
      IF (IREDO .EQ. 0) GO TO 690
C
C     This section computes the predicted values by effectively
C     multiplying the Y array by the pascal triangle matrix.
C     RC is the ratio of new to old values of the coefficient  H*EL(1).
C     when RC differs from 1 by more than 30 percent, or the caller has
C     changed MITER, IWEVAL is set to MITER to force the partials to be
C     updated, if partials are used.  In any case, the partials
C     are updated at least every 20-th step.
C
 200  IF (ABS(RC-1.D0) .GT. 0.3D0) IWEVAL = MITER
      IF (NSTEP .GE. NSTEPJ+20) IWEVAL = MITER
      T = T + H
      DO 210 J1 = 1,NQ
        DO 210 J2 = J1,NQ
          J = (NQ + J1) - J2
          DO 210 I = 1,N
  210       Y(I,J) = Y(I,J) + Y(I,J+1)
C
C     Up to 3 corrector iterations are taken.  A convergence test is
C     made on the r.m.s. norm of each correction, using BND, which
C     is dependent on EPS.  The sum of the corrections is accumulated
C     in the vector ERROR(I).  The Y array is not altered in the
C     corrector loop.  The updated Y vector is stored temporarily in
C     SAVE1.
C
 220  DO 230 I = 1,N
 230    ERROR(I) = 0.D0
      M = 0
      CALL DIFFUN (T, Y, SAVE2, ISTATS)

      IF (MOD(ISTATS,2) .NE. 1) RETURN
      NFE = NFE + 1
      IF (IWEVAL .LE. 0) GO TO 290
C
C     If indicated, the matrix P = I - H*EL(1)*J is reevaluated before
C     starting the corrector iteration.  IWEVAL is set to 0 as an
C     indicator that this has been done.  If MITER = 1 or 2, P is
C     computed and processed in HGPSET.  If MITER = 3, the matrix used
C     is P = I - H*EL(1)*D, where D is a diagonal matrix.
C
      IWEVAL = 0
      RC = 1.D0
      NJE = NJE + 1
      NSTEPJ = NSTEP
      GO TO (250, 240, 260), MITER
 240  NFE = NFE + MW
 250  CON = -H*EL(1)
      CALL HGPSET (DIFFUN, JAC, ISTATS, Y, N0, CON, MITER, IER,
     $      YMAX, ERROR, SAVE1, SAVE2, PW, IPIV)
      IF (MOD(ISTATS,2) .NE. 1) RETURN
      IF (IER .NE. 0) GO TO 420
      GO TO 350
 260  R = EL(1)*.1D0
      DO 270 I = 1,N
 270    PW(I) = Y(I,1) + R*(H*SAVE2(I) - Y(I,2))
      CALL DIFFUN (T, PW, SAVE1, ISTATS)
      IF (MOD(ISTATS,2) .NE. 1) RETURN
      NFE = NFE + 1
      HL0 = H*EL(1)
      DO 280 I = 1,N
        R0 = H*SAVE2(I) - Y(I,2)
        PW(I) = 1.D0
        D = .1D0*R0 - H*(SAVE1(I) - SAVE2(I))
        SAVE1(I) = 0.D0
        IF (ABS(R0) .LT. UROUND*YMAX(I)) GO TO 280
        IF (ABS(D) .EQ. 0.D0) GO TO 420
        PW(I) = .1D0*R0/D
        SAVE1(I) = PW(I)*R0
 280    CONTINUE
      GO TO 370
 290  IF (MITER .NE. 0) GO TO (350, 350, 310), MITER
C
C     In the case of functional iteration, update Y directly from
C     the result of the last DIFFUN call.
C
      D = 0.D0
      DO 300 I = 1,N
        R = H*SAVE2(I) - Y(I,2)
C       ED Apr 1996
C       D = D + ( (R-ERROR(I))/YMAX(I) )**2
        R0 = ABS(R - ERROR(I))
        R1 = MAX(ABS(R), ABS(ERROR(I)))
        IF (R0 .NE. 0.D0 .AND. R1 .NE. 0.D0) THEN
          IF (R0/R1 .LT. UROUND) R0 = 0.D0
        END IF
        D = MAX(D, ABS(R0/YMAX(I)) )
        SAVE1(I) = Y(I,1) + EL(1)*R
 300    ERROR(I) = R
      GO TO 400
C
C     In the case of the chord method, compute the corrector error,
C     F sub (M), and solve the linear system with that as right-hand
C     side and P as coefficient matrix, using the LU decomposition
C     if MITER = 1 or 2.  If MITER = 3, the coefficient H*EL(1)
C     in P is updated.
C
 310  PHL0 = HL0
      HL0 = H*EL(1)
      IF (HL0 .EQ. PHL0) GO TO 330
      R = HL0/PHL0
      DO 320 I = 1,N
        D = 1.D0 - R*(1.D0 - 1.D0/PW(I))
        IF (ABS(D) .EQ. 0.D0) GO TO 440
 320    PW(I) = 1.D0/D
 330  DO 340 I = 1,N
 340    SAVE1(I) = PW(I)*(H*SAVE2(I) - (Y(I,2) + ERROR(I)))
      GO TO 370
 350  DO 360 I = 1,N
 360    SAVE1(I) = H*SAVE2(I) - (Y(I,2) + ERROR(I))
      CALL HGSOL (N0, N, ML, MU, PW, SAVE1, IPIV)
 370  D = 0.D0
      DO 380 I = 1,N
        ERROR(I) = ERROR(I) + SAVE1(I)
C       ED Apr 1996
C       D = D + (SAVE1(I)/YMAX(I))**2
        D = MAX(D, ABS(SAVE1(I)/YMAX(I)) )
        SAVE1(I) = Y(I,1) + EL(1)*ERROR(I)
  380 CONTINUE
C
C     Test for convergence.  If M.GT.0, an estimate of the convergence
C     rate constant is stored in CRATE, and this is used in the test.
C
 400  CONTINUE
C     ED Apr 96
C     IF (M .NE. 0) CRATE = MAX(ZP9*CRATE,D/D1)
      IF (M .NE. 0 .AND. D1 .NE. 0.D0) CRATE = MAX(ZP9*CRATE,D/D1)
      IF ((D*MIN(1.D0,2.D0*CRATE)) .LE. BND) GO TO 450
      D1 = D
      M = M + 1
      IF (M .EQ. 3) GO TO 410
      CALL DIFFUN (T, SAVE1, SAVE2, ISTATS)
      IF (MOD(ISTATS,2) .NE. 1) RETURN
      GO TO 290
C
C     The corrector iteration failed to converge in 3 tries. If partials
C     are involved but are not up to date, they are reevaluated for the
C     next try.  Otherwise the Y array is retracted to its values
C     before prediction, and H is reduced, If possible.  if not, a
C     no-convergence exit is taken.
C
 410  NFE = NFE + 2
      IF (IWEVAL .EQ. -1) GO TO 440
 420  T = TOLD
      RMAX = 2.D0
      DO 430 J1 = 1,NQ
        DO 430 J2 = J1,NQ
          J = (NQ + J1) - J2
          DO 430 I = 1,N
 430        Y(I,J) = Y(I,J) - Y(I,J+1)
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 680
      RH = .25D0
      IREDO = 1
      GO TO 170
 440  IWEVAL = MITER
      GO TO 220
C
C     The corrector has converged.  IWEVAL is set to -1 if partial
C     derivatives were used, to signal that they may need updating on
C     subsequent steps.  The error test is made and control passes to
C     statement 500 if it fails.
C
 450  IF (MITER .NE. 0) IWEVAL = -1
      NFE = NFE + M
      D = 0.D0
      DO 460 I = 1,N
C       ED Apr 1996
C       D = D + (ERROR(I)/YMAX(I))**2
        D = MAX(D, ABS(ERROR(I)/YMAX(I)) )
  460 CONTINUE
      IF (D .GT. E) GO TO 500
C
C     After a successful step, update the Y array.
C     Consider changing H if IDOUB = 1.  Otherwise decrease IDOUB by 1.
C     If IDOUB is then 1 and NQ .LT. MAXDER, then error is saved for
C     use in a possible order increase on the next step.
C     If a change in H is considered, an increase or decrease in order
C     by one is considered also.  A change in H is made only if it is by
C     a factor of at least 1.1.  If not, IDOUB is set to 10 to prevent
C     testing for that many steps.
C
      KFLAG = 0
      IREDO = 0
      NSTEP = NSTEP + 1
      HUSED = H
      NQUSED = NQ
      DO 470 J = 1,L
        DO 470 I = 1,N
 470      Y(I,J) = Y(I,J) + EL(J)*ERROR(I)
      IF (IDOUB .EQ. 1) GO TO 520
      IDOUB = IDOUB - 1
      IF (IDOUB .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    Y(I,LMAX) = ERROR(I)
      GO TO 700
C
C     The error test failed.  KFLAG keeps track of multiple failures.
C     Restore T and the Y array to their previous values, and prepare
C     to try the step again.  Compute the optimum step size for this or
C     one lower order.
C
 500  KFLAG = KFLAG - 1
      T = TOLD
      DO 510 J1 = 1,NQ
        DO 510 J2 = J1,NQ
          J = (NQ + J1) - J2
          DO 510 I = 1,N
 510        Y(I,J) = Y(I,J) - Y(I,J+1)
      RMAX = 2.D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      PR3 = 1.D+20
      GO TO 540
C
C     Regardless of the success or failure of the step, factors
C     PR1, PR2, and PR3 are computed, by which H could be divided
C     at order NQ - 1, order NQ, or order NQ + 1, respectively.
C     In the case of failure, PR3 = 1.E20 to avoid an order increase.
C     The smallest of these is determined and the new order chosen
C     accordingly.  If the order is to be increased, we compute one
C     additional scaled derivative.
C
 520  PR3 = 1.D+20
      IF (L .EQ. LMAX) GO TO 540
      D1 = 0.D0
      DO 530 I = 1,N
C       ED Apr 1996
C       D1 = D1 + ((ERROR(I) - Y(I,LMAX))/YMAX(I))**2
        R0 = ABS(ERROR(I) - Y(I,LMAX))
        R1 = MAX(ABS(ERROR(I)), ABS(Y(I,LMAX)))
        IF (R0 .NE. 0.D0 .AND. R1 .NE. 0.D0) THEN
          IF (R0/R1 .LT. UROUND) R0 = 0.D0
        END IF
        D1 = MAX(D1, ABS(R0/YMAX(I)) )
  530 CONTINUE
C     ED Apr 1996
C     ENQ3 = .5D0/FLOAT(L+1)
      ENQ3 = 1.D0/FLOAT(L+1)
      PR3 = ((D1/EUP)**ENQ3)*1.4D0 + 1.4D-6
  540 CONTINUE
C     ED Apr 1996
C     ENQ2 = .5D0/FLOAT(L)
      ENQ2 = 1.D0/FLOAT(L)
      PR2 = ((D/E)**ENQ2)*1.2D0 + 1.2D-6
      PR1 = 1.D+20
      IF (NQ .EQ. 1) GO TO 560
      D = 0.D0
      DO 550 I = 1,N
C       ED Apr 1996
C       D = D + (Y(I,L)/YMAX(I))**2
        D = MAX(D, ABS(Y(I,L)/YMAX(I)))
  550 CONTINUE
C     ED Apr 1996
C     ENQ1 = .5D0/FLOAT(NQ)
      ENQ1 = 1.D0/FLOAT(NQ)
      PR1 = ((D/EDN)**ENQ1)*1.3D0 + 1.3D-6
C
C     ED Apr 1996
C     Because the computation of D and D1 may involve
C     a large loss of significant digists when ERROR and Y are very
C     close, the estimators PR1, PR2 and PR3 are not very accurate.
C     They do not need to be. Therefore we round RH to 1 digit.
C     This way the algorithm does not propagate
C     this inaccuracy which causes different compilers and different
C     hardware to produce slightly different results.
C
 560  IF (PR2 .LE. PR3) GO TO 570
      IF (PR3 .LT. PR1) GO TO 590
      GO TO 580
 570  IF (PR2 .GT. PR1) GO TO 580
      NEWQ = NQ
      RH = 1.D0/PR2
      J = INT(LOG10(RH))
      I = INT(RH*10.D0**(-J+1)+.5D0)
      RH = 10.D0**(J-1)
      RH = I*RH
      GO TO 620
 580  NEWQ = NQ - 1
      RH = 1.D0/PR1
      J = INT(LOG10(RH))
      I = INT(RH*10.D0**(-J+1)+.5D0)
      RH = 10.D0**(J-1)
      RH = I*RH
      GO TO 620
 590  NEWQ = L
      RH = 1.D0/PR3
      J = INT(LOG10(RH))
      I = INT(RH*10.D0**(-J+1)+.5D0)
      RH = 10.D0**(J-1)
      RH = I*RH
      IF (RH .LT. 1.1D0) GO TO 610
      DO 600 I = 1,N
 600    Y(I,NEWQ+1) = ERROR(I)*EL(L)/FLOAT(L)
      GO TO 630
 610  IDOUB = 10
      GO TO 700
 620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GO TO 610
C
C     If there is a change of order, reset NQ, L, and the coefficients.
C     In any case H is reset according to RH and the Y array is rescaled
C     Then exit from 690 if the step was OK, or redo the step otherwise.
C
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 130
C
C     Control reaches this section if 3 or more failures have occured.
C     It is assumed that the derivatives that have accumulated in the
C     Y array have errors of the wrong order.  Hence the first
C     derivative is recomputed, and the order is set to 1.  Then
C     H is reduced by a factor of 10, and the step is retried.
C     After a total of 7 failures, an exit is taken with KFLAG = -2.
C
 640  IF (KFLAG .EQ. -7) GO TO 670
      RH = .1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      CALL DIFFUN (T, Y, SAVE1, ISTATS)
      IF (MOD(ISTATS,2) .NE. 1) RETURN
      NFE = NFE + 1
      DO 650 I = 1,N
 650    Y(I,2) = H*SAVE1(I)
      IWEVAL = MITER
      IDOUB = 10
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 130
C
C     All returns are made through this section.  H is saved in HOLD
C     to allow the caller to change H on the next step.
C
 660  KFLAG = -1
      GO TO 700
 670  KFLAG = -2
      GO TO 700
 680  KFLAG = -3
      GO TO 700
 690  RMAX = 10.D0
 700  HOLD = H
      JSTART = NQ
      RETURN
      END
