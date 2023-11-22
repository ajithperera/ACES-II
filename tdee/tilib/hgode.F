      SUBROUTINE HGODE (DIFFUN, JAC, NEQN, T0, TSTP, H0, Y0, TOUT, EPS,
     $                  MF, ML, MU, MAXSTP, DT, YP, ISTATS)
C
C     Description of routine.
C
C     Memory manager for Hindmarsh-Gear solver.
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

      INTEGER         I2PRB ( LPI2 )
      COMMON /CI2PRB/ I2PRB

      INTEGER                    LRYFRZ
      EQUIVALENCE ( I2PRB (28) , LRYFRZ )
C
      INCLUDE "comiou.h"
C
C     Declaration of arguments.
C
      INTEGER NEQN, MF, ML, MU, MAXSTP, ISTATS, I
      DOUBLE PRECISION
     $   Y0(NEQN), YP(NEQN), T0, TSTP, H0, TOUT, EPS, DT
      EXTERNAL DIFFUN, JAC
C
C     Declaration of local variables.
C
      INTEGER N, L
C
C     Declaration of local arguments.
C
      INTEGER IPIV(NEQN)
      DOUBLE PRECISION Y(NEQN,13), YMAX(NEQN), ERROR(NEQN),
     $       SAVE1(NEQN), SAVE2(NEQN), PW(NEQN*(2*ML+MU+1))
C
C-----------------------------------------------------------------------
      CALL HGODE0 (N, L, NEQN, MF, ML, MU)
      IF (MOD(ISTATS,2) .EQ. 1) THEN
        IF (COREDY) THEN
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-HGODE  Entering HGODE1"
CSSS                  H0, T0, TOUT,
CSSS     &            EPS, MF, ML, MU, MAXSTP, DT
CSSS      Write(6,*) 
CSSS      Write(6,"(4F10.5)") (Y0(I), I=1,NEQN)
CSSS      Write(6,*) 
CSSS      Write(6,"(4F10.5)") (YP(I), I=1,NEQN)
#endif
           CALL HGODE1 (DIFFUN, JAC, NEQN, T0, TSTP,   H0,  Y0, TOUT,
     $                     EPS, MF,  ML,   MU, MAXSTP, DT,  YP, ISTATS,
     $                      Y, YMAX,ERROR,SAVE1, SAVE2, PW, IPIV)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-HGODE Returning from HGODE1"
#endif

        ELSE
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-HGODE  Entering DIFFUN"
#endif
          CALL DIFFUN(T0,Y0,YP,ISTATS)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-HGODE  Returning DIFFUN"
#endif
          LDOPRP = LMKPRP

          IF (LDOPRP) THEN
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-HGODE  Entering DIFFUN 4DORPP "
#endif
            CALL DIFFUN(T0,Y0,YP,ISTATS)
            LDOPRP = .FALSE.
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "@-HGODE  Returning DIFFUN 4DORPP"
#endif
          END IF
C
C         Analytic Jacobian not implementwed yet. No memory management
C         call required.
C
C         CALL JAC(NEQN,T0,Y0,PD,IPD,ML,MU,ISTATS)
        ENDIF
      ENDIF
      IF (MOD(ISTATS,2) .NE. 1) THEN
        WRITE (IUERR,*) ' @HGODE-F,',' Subroutine exits with error.'
      ENDIF
      RETURN
      END
