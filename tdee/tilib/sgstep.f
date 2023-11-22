      SUBROUTINE  SGSTEP (X,Y,F,ISTATS,NEQN,H,EPS,WT,START,HOLD,K,KOLD,
     $                    CRASH,PHI,P,YP,PSI)
C
C     Description of routine.
C
C     Differential equation solver routines for Adams-Bashforth method.
C     SGSTEP uses a modified divided difference form of
C     the Adams PECE (Predict-Evaluate-Correct-Evaluate) formulas
C     and local extrapolation.  It adjusts the order
C     and step size to control the local error.  Normally each call
C     to SGSTEP advances the solution one step in the direction of TOUT.
C     for reasons of efficiency DE integrates beyond TOUT internally,
C     though never beyond T + 10*(TOUT-T), and calls SGNTRP to
C     interpolate the solution at TOUT. An option is provided to stop
C     the integration at TOUT but it should be used only if it is
C     impossible to continue the integration beyond TOUT .
C
C     Original code: copyright Shampine and Gordon 1975.
C     Modifications: copyright Erik Deumens and QTP, 1990, 1996.
C
      IMPLICIT NONE

C TVG discovered save dependence
      save

      DOUBLE PRECISION
     $     TWOU, FOURU, TEN, TENTH,
     $     ZERO, HALF, SIXTN, HUNDRD, PTNINE, DOS, ONE, FOURTH
      PARAMETER ( TWOU = 1.4D-14, FOURU = 2.8D-14 )
      PARAMETER ( ZERO = 0.D0, HALF = 0.5D0, SIXTN = 16.D0,
     $     HUNDRD = 100.D0,
     $     PTNINE = 0.9D0, DOS = 2.D0, ONE = 1.D0, FOURTH = 0.25D0,
     $     TEN = 10.D0, TENTH = .1D0)
C
C     Arguments
C
      EXTERNAL    F
      LOGICAL START,CRASH
      INTEGER ISTATS, K, KOLD, NEQN
      DOUBLE PRECISION
     $     Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12),
     $     X, H, EPS, HOLD
C
C     Local variables.
C
      LOGICAL PHASE1,NORND
      INTEGER I, L, IQ, J, IP1, KNEW, IFAIL
      INTEGER KP1, KM1, KP2, KM2, NS, NSP1, NSM2, NSP2, IM1, LIMIT1
      DOUBLE PRECISION
     $   ALPHA(12),BETA(12),SIG(13),W(12),V(12),G(13),
     $   GSTR(13),TWO(13)
      DOUBLE PRECISION
     $     XOLD, RHO, HNEW, LIMIT2, REALNS, R, SUM,
     $     ABSH, TAU, ERKM1, ERKM2, ERK, ERKP1, ERR,
     $     TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, REALI, P5EPS,
     $     ROUND
      SAVE ALPHA,BETA,SIG,W,V,G,GSTR,TWO,XOLD
C
C     These machine dependent values are chosen for the Cray.
C     Sun and RS/6000 could have smaller numbers but the
C     difference is not important.
C
      DATA
     $   TWO/2.D0,4.D0,8.D0,16.D0,32.D0,64.D0,128.D0,256.D0,
     $   512.D0,1024.D0,2048.D0,4096.D0,8192.D0/
      DATA
     $   GSTR/5.D-1,8.33D-2,4.17D-2,2.64D-2,1.88D-2,1.43D-2,
     $   1.14D-2,9.36D-3,7.89D-3,6.79D-3,5.92D-3,5.24D-3,
     $   4.68D-3/
C-----------------------------------------------------------------------
C     Check if step size or error tolerance is too small for machine
C     precision.  If first step, initialize PHI array and estimate a
C     starting step size.
C
C     If step size is too small, determine an acceptable one.
C
      CRASH=.TRUE.
      IF ( ABS(H) .LT. FOURU*ABS(X) ) THEN
        H=SIGN(FOURU*ABS(X),H)
        RETURN
      END IF
      P5EPS = HALF*EPS
C
C     If error tolerance is too small, increase it to an acceptable
C     value.
C
      ROUND = ZERO
      DO 10 L=1,NEQN
C       ED Apr 1996
C       ROUND = ROUND + (Y(L)/WT(L))**2
        ROUND = MAX(ROUND, ABS( Y(L)/WT(L) ))
   10 CONTINUE
C     ROUND = TWOU * SQRT(ROUND/FLOAT(NEQN))
      ROUND = TWOU * ROUND
      IF ( P5EPS .LT. ROUND ) THEN
        EPS = DOS*ROUND*(ONE+FOURU)
        RETURN
      END IF
      CRASH=.FALSE.
      G(1) = 1.0D0
      G(2) = 0.5D0
      SIG(1) = 1.0D0
C
      IF ( START ) THEN
C
C       Initialize.  Compute appropriate step size for first step.
C
        CALL F(X,Y,YP,ISTATS)
        IF (MOD(ISTATS,2) .NE. 1) RETURN
        SUM = ZERO
        CALL DCOPY(NEQN,YP,1,PHI(1,1),1)
        CALL DSCAL(NEQN,ZERO,PHI(1,2),1)
        DO 20 L=1,NEQN
C         ED Apr 1996
C         SUM = SUM + (YP(L)/WT(L))**2
          SUM = MAX(SUM, ABS( YP(L)/WT(L) ))
   20   CONTINUE
C       SUM = SQRT(SUM/FLOAT(NEQN))
        ABSH = ABS(H)
        IF (EPS .LT. SIXTN*SUM*H*H) ABSH = FOURTH*SQRT(EPS/SUM)
        H = SIGN(MAX(ABSH,FOURU*ABS(X)),H)
        HOLD = ZERO
        K = 1
        KOLD = 0
        START = .FALSE.
        PHASE1 = .TRUE.
        NORND = .TRUE.
        IF ( P5EPS .LE. HUNDRD*ROUND ) THEN
          NORND = .FALSE.
          CALL DSCAL(NEQN,ZERO,PHI(1,15),1)
        END IF
      END IF
      IFAIL = 0
C
C     Compute coefficients of formulas for this step.  Avoid computing
C     those quantities not changed when step size is not changed.
C
  100 CONTINUE
      KP1 = K+1
      KP2 = K+2
      KM1 = K-1
      KM2 = K-2
C
C     NS is the number of steps taken with size h, including the current
C     one.  When K.LT.NS, no coefficients change.
C
      IF (H.NE.HOLD) NS = 0
      IF (NS.LE.KOLD) NS = NS+1
      NSP1 = NS+1
C
      IF ( K .GE. NS ) THEN
C
C       Compute those components of ALPHA(*),BETA(*),PSI(*),SIG(*) which
C       are changed.
C
        BETA(NS) = ONE
        REALNS = NS
        ALPHA(NS) = ONE/REALNS
        TEMP1 = H*REALNS
        SIG(NSP1) = ONE
        IF ( K .GE. NSP1 ) THEN
          DO 105 I = NSP1,K
            IM1 = I-1
            TEMP2 = PSI(IM1)
            PSI(IM1) = TEMP1
            BETA(I) = BETA(IM1)*PSI(IM1)/TEMP2
            TEMP1 = TEMP2+H
            ALPHA(I) = H/TEMP1
            REALI = I
            SIG(I+1) = REALI*ALPHA(I)*SIG(I)
  105     CONTINUE
        END IF
        PSI(K) = TEMP1
C
C       Compute coefficients g(*)
C
C       Initialize V(*) and set W(*).  G(2) is set before.
C
        IF ( NS .LE. 1 ) THEN
          DO 115 IQ = 1,K
            TEMP3 = IQ*(IQ+1)
            V(IQ) = ONE/TEMP3
            W(IQ) = V(IQ)
  115     CONTINUE
        ELSE
C
C         If order was raised, update diagonal part of V(*).
C
          IF ( K .GT. KOLD) THEN
            TEMP4 = K*KP1
            V(K) = ONE/TEMP4
            NSM2 = NS-2
            IF  ( NSM2 .GE. 1 ) THEN
              DO 125 J = 1,NSM2
                I = K-J
                V(I) = V(I)-ALPHA(J+1)*V(I+1)
  125         CONTINUE
            END IF
          END IF
C
C         Update V(*) and set W(*).
C
          LIMIT1 = KP1-NS
          TEMP5 = ALPHA(NS)
          DO 135 IQ = 1,LIMIT1
            V(IQ) = V(IQ)-TEMP5*V(IQ+1)
            W(IQ) = V(IQ)
  135     CONTINUE
          G(NSP1) = W(1)
        END IF
C
C       Compute the G(*) in the work vector W(*).
C
        NSP2 = NS+2
        IF ( KP1 .GE. NSP2 ) THEN
          DO 150 I = NSP2,KP1
            LIMIT2 = KP2-I
            TEMP6 = ALPHA(I-1)
            DO 145 IQ = 1,INT(LIMIT2)
              W(IQ) = W(IQ)-TEMP6*W(IQ+1)
  145       CONTINUE
            G(I) = W(1)
  150     CONTINUE
        END IF
      END IF
C
C     Predict a solution P(*), evaluate derivatives using predicted
C     solution, estimate local error at order k and errors at orders K,
C     K-1, K-2 as if constant step size were used.
C
C     Change PHI to PHI star.
C
      IF ( K .GE. NSP1 ) THEN
        DO 210 I = NSP1,K
           CALL DSCAL(NEQN,BETA(I),PHI(1,I),1)
  210   CONTINUE
      END IF
C
C     Predict solution and differences.
C
        CALL DCOPY(NEQN,PHI(1,KP1),1,PHI(1,KP2),1)
        CALL DSCAL(NEQN,ZERO,PHI(1,KP1),1)
        CALL DSCAL(NEQN,ZERO,P,1)
      DO 230 J = 1,K
        I = KP1-J
        IP1 = I+1
        TEMP2 = G(I)
        CALL DAXPY(NEQN,G(I),PHI(1,I),1,P,1)
        CALL DAXPY(NEQN,ONE,PHI(1,IP1),1,PHI(1,I),1)
  230 CONTINUE
      IF ( .NOT. NORND ) THEN
        DO 235 L = 1,NEQN
          TAU = H*P(L)-PHI(L,15)
          P(L) = Y(L)+TAU
          PHI(L,16) = (P(L)-Y(L))-TAU
  235   CONTINUE
      ELSE
        CALL DSCAL(NEQN,H,P,1)
        CALL DAXPY(NEQN,ONE,Y,1,P,1)
      END IF
      XOLD = X
      X = X+H
      ABSH = ABS(H)

C      print*,"Current order =",k,"relerr =",eps,"step =",H

      CALL F(X,P,YP,ISTATS)
      IF (MOD(ISTATS,2) .NE. 1) RETURN
C
C     Estimate errors at orders K,K-1,K-2.
C
      ERKM2 = ZERO
      ERKM1 = ZERO
      ERK = ZERO
      DO 265 L = 1,NEQN
        TEMP3 = ONE/WT(L)
        TEMP4 = YP(L)-PHI(L,1)
        IF ( KM2 .GT. 0 ) THEN
C         ED Apr 1996
C         ERKM2 = ERKM2 + ((PHI(L,KM1)+TEMP4)*TEMP3)**2
          ERKM2 = MAX(ERKM2, ABS( (PHI(L,KM1)+TEMP4)*TEMP3 ))
        END IF
        IF ( KM2 .GE. 0 ) THEN
C         ED Apr 1996
C         ERKM1 = ERKM1 + ((PHI(L,K)+TEMP4)*TEMP3)**2
          ERKM1 = MAX(ERKM1, ABS( (PHI(L,K)+TEMP4)*TEMP3 ))
        END IF
C       ED Apr 1996
C       ERK = ERK + (TEMP4*TEMP3)**2
        ERK = MAX(ERK, ABS( TEMP4*TEMP3 ))
  265 CONTINUE
C
      IF ( KM2 .GT. 0 ) THEN
C       ED Apr 1996
C       ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*SQRT(ERKM2/FLOAT(NEQN))
        ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*ERKM2
      END IF
      IF ( KM2 .GE. 0 ) THEN
C       ED Apr 1996
C       ERKM1 = ABSH*SIG(K)*GSTR(KM1)*SQRT(ERKM1/FLOAT(NEQN))
        ERKM1 = ABSH*SIG(K)*GSTR(KM1)*ERKM1
      END IF
C     ED Apr 1996
C     TEMP5 = ABSH * SQRT(ERK/FLOAT(NEQN))
      TEMP5 = ABSH * ERK
      ERR = TEMP5 * (G(K)-G(KP1))
      ERK = TEMP5 * SIG(KP1)*GSTR(K)
      KNEW = K
C
C     Test if order should be lowered.
C
      IF ( KM2 .GT. 0 ) THEN
        IF ( MAX(ERKM1,ERKM2) .LE. ERK ) KNEW = KM1
      ELSE IF ( KM2 .EQ. 0 ) THEN
        IF ( ERKM1 .LE. HALF*ERK ) KNEW = KM1
      END IF
C
C     Test if step successful.
C
      IF ( ERR .GT. EPS ) THEN
C
C       The step is unsuccessful.  Restore  X, PHI(*,*), PSI(*) .
C       if third consecutive failure, set order to one.  if step fails
C       more than three times, consider an optimal step size.  Double
C       error tolerance and return if estimated step size is too small
C       for machine precision.
C
C       Restore X, PHI(*,*) and PSI(*).
C
        PHASE1 = .FALSE.
        X = XOLD
        DO 310 I = 1,K
          TEMP1 = ONE/BETA(I)
          IP1 = I+1
          CALL DAXPY(NEQN,-ONE,PHI(1,IP1),1,PHI(1,I),1)
          CALL DSCAL(NEQN,TEMP1,PHI(1,I),1)
  310   CONTINUE
        IF ( K .GE. 2 ) THEN
          DO 315 I = 2,K
            PSI(I-1) = PSI(I)-H
  315     CONTINUE
        END IF
C
C       On third failure, set order to one.  thereafter, use optimal
C       step size.
C
        IFAIL = IFAIL+1
        TEMP2 = HALF
C
        IF ( IFAIL-3 .GT. 0) THEN
          IF (P5EPS.LT.FOURTH*ERK) TEMP2 = SQRT(P5EPS/ERK)
        END IF
        IF ( IFAIL-3 .GE. 0) THEN
          KNEW = 1
        END IF
        H = TEMP2*H
        K = KNEW
        IF ( ABS(H) .LT. FOURU*ABS(X) ) THEN
          CRASH = .TRUE.
          H = SIGN(FOURU*ABS(X),H)
          EPS = EPS+EPS
          RETURN
        END IF
        GO TO 100
      END IF
C
C     The step is successful.  correct the predicted solution, evaluate
C     the derivatives using the corrected solution and update the
C     differences.  determine best order and step size for next step.
C
      KOLD = K
      HOLD = H
      TEMP1 = H*G(KP1)
C
C     Correct and evaluate.
C
      IF (.NOT. NORND) THEN
        DO 405 L = 1,NEQN
          RHO = TEMP1*(YP(L)-PHI(L,1))-PHI(L,16)
          Y(L) = P(L)+RHO
          PHI(L,15) = (Y(L)-P(L))-RHO
  405   CONTINUE
      ELSE
        CALL DCOPY(NEQN,YP,1,Y,1)
        CALL DAXPY(NEQN,-ONE,PHI(1,1),1,Y,1)
        CALL DSCAL(NEQN,TEMP1,Y,1)
        CALL DAXPY(NEQN,ONE,P,1,Y,1)
      END IF
C
      CALL F(X,Y,YP,ISTATS)
      IF (MOD(ISTATS,2) .NE. 1) RETURN
C
C     Update differences for next step.
C
      CALL DCOPY(NEQN,PHI(1,1),1,PHI(1,KP1),1)
      CALL DSCAL(NEQN,-ONE,PHI(1,KP1),1)
      CALL DAXPY(NEQN,ONE,YP,1,PHI(1,KP1),1)
      CALL DSCAL(NEQN,-ONE,PHI(1,KP2),1)
      CALL DAXPY(NEQN,ONE,PHI(1,KP1),1,PHI(1,KP2),1)
      DO 435 I = 1,K
         CALL DAXPY(NEQN,ONE,PHI(1,KP1),1,PHI(1,I),1)
  435 CONTINUE
C
C     Estimate error at order k+1 unless:
C     in first phase when always raise order,
C     already decided to lower order,
C     step size not constant so estimate unreliable
C
      ERKP1 = ZERO
      IF (KNEW.EQ.KM1 .OR. K.EQ.12) PHASE1 = .FALSE.
C
      IF ( .NOT. (PHASE1 .OR. KNEW.EQ.KM1 .OR. KP1.GT.NS) ) THEN
        DO 440 L = 1,NEQN
C         ED Apr 1996
C         ERKP1 = ERKP1+(PHI(L,KP2)/WT(L))**2
          ERKP1 = MAX(ERKP1, ABS( PHI(L,KP2)/WT(L) ))
  440   CONTINUE
C       ERKP1 = ABSH * GSTR(KP1) * SQRT(ERKP1/FLOAT(NEQN))
        ERKP1 = ABSH * GSTR(KP1) * ERKP1
C
C       Using estimated error at order k+1, determine appropriate order
C       for next step.
C
        IF ( (K.LE.1 .AND. ERKP1.LT.HALF*ERK)
     $     .OR.
     $     ( K.GT.1 .AND. ERKM1.GT.MIN(ERK,ERKP1)
     $       .AND. (ERKP1.LT.ERK .AND. K.NE.12)) ) THEN
C
C         Here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would
C         have been lowered in block 2.  Thus order is to be raised.
C
C         Raise order.
C
          K = KP1
          ERK = ERKP1
        ELSE IF ( K.GT.1 .AND. ERKM1.LE.MIN(ERK,ERKP1) ) THEN
C
C         Lower order.
C
          K = KM1
          ERK = ERKM1
        END IF
      ELSE IF ( PHASE1 ) THEN
C
C       Raise order.
C
        K = KP1
        ERK = ERKP1
      ELSE IF ( KNEW .EQ. KM1 ) THEN
C
C       Lower order.
C
        K = KM1
        ERK = ERKM1
      END IF
C
C     With new order determine appropriate step size for next step.
C
      HNEW = H+H
C
      IF ( .NOT.PHASE1 .AND. P5EPS.LT.ERK*TWO(K+1) ) THEN
        HNEW = H
        IF ( P5EPS .LT. ERK ) THEN
          TEMP2 = K+1
          R = (P5EPS/ERK)**(ONE/TEMP2)
C
C         ED Apr 1996 >>>
C         Because the computation of P5EPS and ERK and TEMP2 may involve
C         a large loss of significant digits when YP and PHI are very
C         close, the estimator R is not very accurate.
C         It does not need to be. Therefore we round them up to 1 digit
C         before using them. This way the algorithm does not propagate
C         this inaccuracy which causes different compilers and different
C         hardware to produce slightly different results.
C
          J = INT(LOG10(R))
          I = INT(R*TEN**(-J+1)+HALF)
          R = TENTH**(J-1)
          R = R/I
C         ED Apr 1996 <<<
          HNEW = ABSH*MAX(HALF,MIN(PTNINE,R))
          HNEW = SIGN(MAX(HNEW,FOURU*ABS(X)),H)
        END IF
      END IF
      H = HNEW
      RETURN
      END
