      SUBROUTINE HGODE0 (N,L,NEQN,MF,ML,MU)
C
C
C     Description of routine.
C     Calculates scratch block lengths in words for HGODE.
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
      INTEGER N, L, NEQN, MF, ML, MU
      INTEGER I
C----------------------------------------------------------------------
      N = 7
      NAMBLK(1) = 'Y'
      IF ( MF/10 .EQ. 2 ) THEN
        LENIND(1) = NWRDBL * NEQN*6
      ELSE
        LENIND(1) = NWRDBL * NEQN*13
      ENDIF
      NAMBLK(2) = 'YMAX'
      LENIND(2) = NWRDBL * NEQN
      NAMBLK(3) = 'ERROR'
      LENIND(3) = NWRDBL * NEQN
      NAMBLK(4) = 'SAVE1'
      LENIND(4) = NWRDBL * NEQN
      NAMBLK(5) = 'SAVE2'
      LENIND(5) = NWRDBL * NEQN
C
C     Space for real PW and integer IPIV
C
      NAMBLK(6) = 'PW'
      NAMBLK(7) = 'IPIV'
      IF ( MOD(MF, 10) .EQ. 0 ) THEN
        LENIND(6) = NWRDBL * 1
        LENIND(7) = 1
      ELSE IF ( MOD(MF, 10) .EQ. 1 .OR. MOD(MF, 10) .EQ. 2 ) THEN
        LENIND(6) = NWRDBL * NEQN*(2*ML+MU+1)
        LENIND(7) = NEQN
      ELSE IF ( MOD(MF, 10) .EQ. 3 ) THEN
        LENIND(6) = NWRDBL * NEQN
        LENIND(7) = 1
      ENDIF
      L = 0
      DO 100 I=1,N
        L = L + LENIND(I)
  100 CONTINUE
      RETURN
      END
