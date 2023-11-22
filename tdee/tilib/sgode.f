      SUBROUTINE SGODE (F, NEQN, Y, YP, T, TOUT, RELERR, ABSERR,
     $   MAXNUM, DT, ISTATS)
C
C
C     Description of routine.
C
C     Memory manager for Shampine-Gordon solvecr.
C
C     Original code: copyright Erik Deumens and QTP, 1992.
C     Modifications: copyright Erik Deumens and QTP, 1996.
C
      IMPLICIT NONE
C
C     Declaration of global variables.
C
      INCLUDE "lengths.h"
      INCLUDE "comcor.h"
C
      LOGICAL         L1PRB ( LPL1 )
      COMMON /CL1PRB/ L1PRB
      LOGICAL                    LMKPRP
      EQUIVALENCE ( L1PRB (34) , LMKPRP )
C
      LOGICAL         L2PRB ( LPL2 )
      COMMON /CL2PRB/ L2PRB
      LOGICAL                    LDOPRP
      EQUIVALENCE ( L2PRB (12) , LDOPRP )
C
      INCLUDE "comiou.h"
C
C     Declaration of arguments.
C
      INTEGER NEQN, ISTATS, MAXNUM
      DOUBLE PRECISION 
     $  Y(NEQN), YP(NEQN), T, TOUT, DT, RELERR, ABSERR
      EXTERNAL F
C
C     Declaration of local variables.
C
      INTEGER  N, L
C
C     Declaration of arguments.
C
      DOUBLE PRECISION YY(NEQN), WT(NEQN), P(NEQN),
     $                 PHI(NEQN,16), YPOUT(NEQN)
C
C---------------------------------------------------------------------
      CALL DZERO(YY,NEQN); CALL DZERO(WT,NEQN); CALL DZERO(P,NEQN)
      CALL DZERO(PHI,NEQN*16); CALL DZERO(YPOUT,NEQN)
      CALL SGODE0 (N, L, NEQN)
      IF (MOD(ISTATS,2) .EQ. 1) THEN
        IF (COREDY) THEN
            CALL SGODE1(F, NEQN, Y, YP, T, TOUT, RELERR, ABSERR,
     $                  MAXNUM, DT, ISTATS, YY, WT, PHI, P, YPOUT)
        ELSE
          CALL F(T,Y,YP,ISTATS)
          LDOPRP = LMKPRP
          IF (LDOPRP) THEN
            CALL F(T,Y,YP,ISTATS)
            LDOPRP = .FALSE.
          END IF
        END IF
      END IF
      IF (MOD(ISTATS,2) .NE. 1) THEN
        WRITE (IUERR,*) ' @SGODE-F,',' Subroutine exits with error.'
      END IF
      RETURN
      END
