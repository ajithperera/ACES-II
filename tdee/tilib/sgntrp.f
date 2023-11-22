      SUBROUTINE SGNTRP (X,Y,XOUT,YOUT,YPOUT,NEQN,KOLD,PHI,PSI)
C
C
C     Description of routine.
C
C     Differential equation solver routines for Adams-Bashforth method.
C
C     The methods in subroutine SGSTEP approximate the solution near X
C     by a polynomial.  Subroutine SGNTRP approximates the solution at
C     XOUT by evaluating the polynomial there.  Information defining
C     this polynomial is passed from SGSTEP so SGNTRP cannot be used
C     alone.
C
C     This code is completely explained and documented in the text,
C     Computer Solution of Ordinary Differential Equations:  The Initial
C     Value Problem, by L. F. Shampine and M. K. Gordon (Freeman, 1975)
C
C     Input to SGNTRP
C
C     The user provides storage in the calling program for the arrays in
C     the call list and defines
C
C     XOUT -- point at which solution is desired.
C
C     The remaining parameters are defined in SGSTEP and passed to
C     SGNTRP from that subroutine
C
C     Output from  SGNTRP
C
C     YOUT(*) -- solution at  XOUT
C
C     YPOUT(*) -- derivative of solution at  XOUT
C
C     The remaining parameters are returned unaltered from their input
C     values.  Integration with SGSTEP may be continued.
C
C     Original code: copyright Shampine and Gordon 1975.
C     Modifications: copyright Erik Deumens and QTP, 1990, 1996.
C
      IMPLICIT NONE
      INTEGER KOLD, NEQN, I, J, JM1, KIP1, KI, LIMIT1
      DOUBLE PRECISION
     $   X, Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),PSI(12),
     $   G(13),W(13),RHO(13),
     $   ONE, ZERO, PSIJM1, GAMMA, ETA,
     $   TEMP1, TERM, XOUT, HI
      SAVE G, W, RHO, ONE, ZERO
      DATA G(1)/1.0D0/,RHO(1)/1.0D0/
      DATA ONE/1.D0/,ZERO/0.D0/
C-----------------------------------------------------------------------
      HI = XOUT - X
      KI = KOLD + 1
      KIP1 = KI + 1
C
C     Initialize W(*) for computing G(*)
C
      DO 5 I = 1,KI
        TEMP1 = I
        W(I) = ONE/TEMP1
    5 CONTINUE
      TERM = ZERO
C
C     Compute G(*)
C
      DO 15 J = 2,KI
        JM1 = J - 1
        PSIJM1 = PSI(JM1)
        GAMMA = (HI + TERM)/PSIJM1
        ETA = HI/PSIJM1
        LIMIT1 = KIP1 - J
        DO 10 I = 1,LIMIT1
          W(I) = GAMMA*W(I) - ETA*W(I+1)
   10   CONTINUE
        G(J) = W(1)
        RHO(J) = GAMMA*RHO(JM1)
        TERM = PSIJM1
   15 CONTINUE
C
C     Interpolate
C
                   CALL DSCAL
     $ (NEQN,ZERO,YOUT,1)
                   CALL DSCAL
     $ (NEQN,ZERO,YPOUT,1)
      DO 30 J = 1,KI
        I = KIP1 - J
                   CALL DAXPY
     $   (NEQN,G(I),PHI(1,I),1,YOUT,1)
                   CALL DAXPY
     $   (NEQN,RHO(I),PHI(1,I),1,YPOUT,1)
   30 CONTINUE
                   CALL DSCAL
     $ (NEQN,HI,YOUT,1)
                   CALL DAXPY
     $ (NEQN,ONE,Y,1,YOUT,1)
      RETURN
      END
