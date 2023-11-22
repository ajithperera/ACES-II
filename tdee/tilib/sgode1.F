      SUBROUTINE SGODE1 (F, NEQN, Y, YP, T, TOUT, RELERR, ABSERR,
     $   MAXSTP, DT, ISTATS, YY, WT, PHI, P, YPOUT)
C
C
C     Description of routine.
C
C     Differential equation solver routines for Adams-Bashforth method.
C     Keeps temporary arrays for continued integration and calls WRTRST.
C
C     Original code: copyright Erik Deumens and QTP, 1992.
C     Modifications: copyright Erik Deumens and QTP, 1996.
C
      IMPLICIT NONE
C
C     Declaration of global variables.
C
      INCLUDE "lengths.h"
C
      LOGICAL         L1PRB ( LPL1 )
      COMMON /CL1PRB/ L1PRB

      LOGICAL                    LMKPRP
      EQUIVALENCE ( L1PRB (34) , LMKPRP )
C
      LOGICAL         L2PRB ( LPL2 )
      COMMON /CL2PRB/ L2PRB

      LOGICAL                    LBKSTP
      EQUIVALENCE ( L2PRB (11) , LBKSTP )

      LOGICAL                    LDOPRP
      EQUIVALENCE ( L2PRB (12) , LDOPRP )
C
      INTEGER         I2PRB ( LPI2 )
      COMMON /CI2PRB/ I2PRB

      INTEGER                    KDYNFU
      EQUIVALENCE ( I2PRB (61) , KDYNFU )

      INTEGER                     KODFN0
      EQUIVALENCE ( I2PRB (113) , KODFN0 )

      INTEGER                     KODSTP
      EQUIVALENCE ( I2PRB (135) , KODSTP )
C
      DOUBLE PRECISION R1PRB ( LPR1 )
      COMMON /CR1PRB/  R1PRB

      DOUBLE PRECISION           RSMALL
      EQUIVALENCE ( R1PRB ( 1) , RSMALL )
C
      DOUBLE PRECISION R2PRB ( LPR2 )
      COMMON /CR2PRB/  R2PRB

      DOUBLE PRECISION           RMNSTP
      EQUIVALENCE ( R2PRB (48) , RMNSTP )

      DOUBLE PRECISION           RMXSTP
      EQUIVALENCE ( R2PRB (49) , RMXSTP )

      DOUBLE PRECISION           RSTPSM
      EQUIVALENCE ( R2PRB (50) , RSTPSM )

      DOUBLE PRECISION           RSTPSQ
      EQUIVALENCE ( R2PRB (51) , RSTPSQ )
C
      INCLUDE "comiou.h"
C
C     Declaration of arguments.
C
      INTEGER NEQN, ISTATS, MAXSTP
      EXTERNAL F
      DOUBLE PRECISION Y(NEQN), YP(NEQN), YY(NEQN), WT(NEQN),
     $                 P(NEQN), YPOUT(NEQN), PHI(NEQN,16)
      DOUBLE PRECISION T, TOUT, DT, RELERR, ABSERR
C
C     Declaration of local variables.
C
      INTEGER I, IFLAG
      DOUBLE PRECISION TSTP, T0
C-----------------------------------------------------------------------
      IF (KODSTP .EQ. 0) KODFN0 = KDYNFU
      IFLAG = 1
      T0 = T
      TSTP = T
      DO 100 I=1,MAXSTP
        TSTP = TSTP + SIGN(DT,TOUT-T0)
C
C       Do not overshoot end point.
C
        IF ((TSTP - TOUT) * (TSTP - T0) .GT. 0.D0) TSTP = TOUT
        CALL SGDRVR (F, NEQN, Y, YP, T, TSTP, RELERR, ABSERR,
     $               MAXSTP, IFLAG,
     $               KODSTP, RMNSTP, RMXSTP, RSTPSM, RSTPSQ, ISTATS,
     $               YY, WT, PHI, P, YPOUT)
        IF (MOD(ISTATS,2) .EQ. 0) THEN
C
C         Integrator failed, try afresh from last restart file.
C
          RETURN
        ELSE
C
C         All is well, call the function on trajectory with properties
C         and write restart file.
C
          IF (IFLAG .EQ. 2) THEN
            LDOPRP = LMKPRP
            CALL F (T,Y,YP,ISTATS)
            IF (LDOPRP) LDOPRP = .FALSE.
#ifdef _DEBUG_LVL0
            Write(6,*)
            Write(6,*) "@-SGODE Restart write"
#endif
            CALL RSTMGR ('Write', ISTATS)
            LBKSTP = .FALSE.
C
C           Exit if the end is reached.
C
            IF ((T-TOUT) * (TOUT-T0) .GE. -RSMALL) RETURN
          ELSE
C
C           Exit on error.
C
            WRITE (IUERR,*) ' @SGODE-F,',
     $         ' IFLAG indicates an error:'
            WRITE (IUERR,*) '           ',' IFLAG = ',IFLAG,'.'
            ISTATS = 0
            RETURN
          ENDIF
        ENDIF
  100 CONTINUE
      WRITE (IUERR,*) ' @SGODE-F,',
     $            ' DE equation solver needed too many steps.'
      WRITE (IUERR,*) '           ',' MAXSTP = ',MAXSTP,'.'
      ISTATS = 0
      RETURN
      END
