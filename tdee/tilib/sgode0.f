      SUBROUTINE SGODE0 (N,L,NEQN)
C
C
C     Description of routine.
C
C     Differential equation solver routines for Adams-Bashforth method.
C     Calculates scratch block lengths in words for SGINT.
C
C     Original code: copyright Erik Deumens and QTP, 1992.
C     Modifications: copyright Erik Deumens and QTP, 1996.
C
      IMPLICIT NONE
C
C     Declaration of global variables.
C
      INCLUDE "comcor.h"
C
C     Declaration of arguments.
C
      INTEGER N, L, NEQN
      INTEGER I
C----------------------------------------------------------------------
      N = 5
      NAMBLK(1) = 'YY'
      LENIND(1) = NWRDBL * NEQN
      NAMBLK(2) = 'WT'
      LENIND(2) = NWRDBL * NEQN
      NAMBLK(3) = 'PHI'
      LENIND(3) = NWRDBL * NEQN*16
      NAMBLK(4) = 'P'
      LENIND(4) = NWRDBL * NEQN
      NAMBLK(5) = 'YPOUT'
      LENIND(5) = NWRDBL * NEQN
      L = 0
      DO 100 I=1,N
        L = L + LENIND(I)
  100 CONTINUE
      RETURN
      END
