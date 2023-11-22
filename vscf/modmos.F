
C This routine reads records EVCSYMAC and EVCSYMBC to see if the AO
C basis MO coefficients transform as irreducible representations of
C the computational point group.

      SUBROUTINE MODMOS(EVEC,SCR,EVAL,NBASX,NBAS,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#include "maxbasfn.par"
      INTEGER NBASX,NBAS,IUHF
      CHARACTER*8 IRRSYM(MAXBASFN)
      DIMENSION EVEC(NBASX,NBAS,IUHF+1)
      DIMENSION SCR(NBASX,NBAS)
      DIMENSION EVAL(NBAS,IUHF+1)

      CHARACTER*4 PTGRP
      INTEGER NPAIR,PAIR1,PAIR2
      LOGICAL MOSYMOK,COUNTED
      DIMENSION MOSYMOK(MAXBASFN,2),COUNTED(MAXBASFN,2)
      DIMENSION NPAIR(2),PAIR1(MAXBASFN,2),PAIR2(MAXBASFN,2)
      DIMENSION CNEW(MAXBASFN)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      DATA DGETOL /0.00001D+00/

      CALL GETCREC(20,'JOBARC','COMPPTGP',4,PTGRP)
C
C     Check if MOs transform as irreducible representations of the
C     computational point group.
C
      DO 20 ISPIN=1,IUHF+1
C
      IF(ISPIN.EQ.1)THEN
       CALL GETREC(20,'JOBARC','EVCSYMAC',IINTFP*NBAS,IRRSYM)
      ELSE
       CALL GETREC(20,'JOBARC','EVCSYMBC',IINTFP*NBAS,IRRSYM)
      ENDIF
C
      DO 10 IMO=1,NBAS
C
      MOSYMOK(IMO,ISPIN) = .TRUE.
C
      IF( IRRSYM(IMO)(1:4) .EQ. '    ' .OR.
     &    IRRSYM(IMO)(1:4) .EQ. '   g' .OR.
     &    IRRSYM(IMO)(1:4) .EQ. '   u'        )THEN
       WRITE(6,1000) ISPIN,IMO
       MOSYMOK(IMO,ISPIN) = .FALSE.
      ENDIF
C
C     Special Cs group code.
C
      IF(PTGRP.EQ.'C s ' .OR. PTGRP.EQ.'C1h ')THEN
       IF(IRRSYM(IMO)(1:4) .EQ. ' A  ')THEN
        WRITE(6,1000) ISPIN,IMO
        WRITE(6,'(A4)') IRRSYM(IMO)(1:4)
        MOSYMOK(IMO,ISPIN) = .FALSE.
       ENDIF
      ENDIF
C
C     Special C2v group code.
C
      IF(PTGRP.EQ.'C2v ')THEN
       IF(IRRSYM(IMO)(1:4) .NE. '  A1' .AND.
     &    IRRSYM(IMO)(1:4) .NE. '  A2' .AND.
     &    IRRSYM(IMO)(1:4) .NE. '  B1' .AND.
     &    IRRSYM(IMO)(1:4) .NE. '  B2')THEN
        WRITE(6,1000) ISPIN,IMO
        WRITE(6,'(A4)') IRRSYM(IMO)(1:4)
        MOSYMOK(IMO,ISPIN) = .FALSE.
       ENDIF
      ENDIF
C
   10 CONTINUE
   20 CONTINUE
C
C     We have array telling us which MOs have correct symmetry. Loop
C     through set and see if any of the ones with incorrect symmetry
C     are nondegenerate. If so, we have problems.
C
      DO 50 ISPIN=1,IUHF+1
      DO 40   IMO=1,NBAS
C
      NDEG = 1
      IF(.NOT. MOSYMOK(IMO,ISPIN))THEN
       DO 30  JMO=1,NBAS
        IF(JMO.NE.IMO)THEN
         IF( DABS(EVAL(IMO,ISPIN) - EVAL(JMO,ISPIN)) .LT. DGETOL)THEN
          NDEG = NDEG + 1
         ENDIF
        ENDIF
   30  CONTINUE
       WRITE(6,1010) IMO,ISPIN,NDEG
C
       IF(NDEG.EQ.1)THEN
        WRITE(6,*) ' @MODMOS-F, Nondegenerate orbital problem. '
        CALL ERREX
       ENDIF
C
       IF(NDEG.GT.2)THEN
        WRITE(6,*) ' @MODMOS-F, Degeneracy too high. '
        CALL ERREX
       ENDIF
C
      ENDIF
C
   40 CONTINUE
   50 CONTINUE
C
C     Now try to match the degenerate pairs.
C
      DO 80 ISPIN=1,IUHF+1
      DO 70   IMO=1,NBAS
       COUNTED(IMO,ISPIN) = .FALSE.
   70 CONTINUE
   80 CONTINUE
C
      DO 110 ISPIN=1,IUHF+1
      IPAIR = 0
      DO 100   IMO=1,NBAS
C
      IF(.NOT. MOSYMOK(IMO,ISPIN) .AND. .NOT. COUNTED(IMO,ISPIN) )THEN
       DO 90  JMO=IMO+1,NBAS
        IF( DABS(EVAL(IMO,ISPIN) - EVAL(JMO,ISPIN)) .LT. DGETOL)THEN
         IPAIR = IPAIR + 1
         PAIR1(IPAIR,ISPIN) = IMO
         PAIR2(IPAIR,ISPIN) = JMO
        ENDIF
   90  CONTINUE
      ENDIF
  100 CONTINUE
      NPAIR(ISPIN) = IPAIR
  110 CONTINUE
C
C     Write out the pair information.
C
      DO 150 ISPIN=1,IUHF+1
      DO 140 IPAIR=1,NPAIR(ISPIN)
C
      IMO = PAIR1(IPAIR,ISPIN)
      JMO = PAIR2(IPAIR,ISPIN)
      WRITE(6,1020) ISPIN,IPAIR,IMO,JMO,EVAL(IMO,ISPIN),EVAL(JMO,ISPIN)
  140 CONTINUE
  150 CONTINUE
C
C     Now try to change the orbitals.
C
      DO 170 ISPIN=1,IUHF+1
C
c YAU : old
c     CALL ICOPY(IINTFP*NBASX*NBAS,EVEC(1,1,ISPIN),1,SCR,1)
c YAU : new
      CALL DCOPY(NBASX*NBAS,EVEC(1,1,ISPIN),1,SCR,1)
c YAU : end
C
      DO 160 IPAIR=1,NPAIR(ISPIN)
C
      IMO = PAIR1(IPAIR,ISPIN)
      JMO = PAIR2(IPAIR,ISPIN)
C
      CALL VADD(CNEW,SCR(1,IMO),SCR(1,JMO),NBASX, 1.0D+00)
      CALL SSCAL(NBASX,DSQRT(2.0D+00),CNEW,1)
c YAU : old
c     CALL ICOPY(IINTFP*NBASX,CNEW,1,EVEC(1,IMO,ISPIN),1)
c YAU : new
      CALL DCOPY(NBASX,CNEW,1,EVEC(1,IMO,ISPIN),1)
c YAU : end
C
      CALL VADD(CNEW,SCR(1,IMO),SCR(1,JMO),NBASX,-1.0D+00)
      CALL SSCAL(NBASX,DSQRT(2.0D+00),CNEW,1)
c YAU : old
c     CALL ICOPY(IINTFP*NBASX,CNEW,1,EVEC(1,JMO,ISPIN),1)
c YAU : new
      CALL DCOPY(NBASX,CNEW,1,EVEC(1,JMO,ISPIN),1)
c YAU : end
C
  160 CONTINUE
C
      IF(ISPIN.EQ.1)THEN
       CALL PUTREC(20,'JOBARC','EVECAO_A',NBASX*NBAS*IINTFP,
     &             EVEC(1,1,ISPIN))
      ELSE
       CALL PUTREC(20,'JOBARC','EVECAO_B',NBASX*NBAS*IINTFP,
     &             EVEC(1,1,ISPIN))
      ENDIF
  170 CONTINUE
C
      RETURN
 1000 FORMAT(' @MODMOS-I, Spin ',I3,' Orbital ',I5,
     &       ' is not symmetry-adapted. ')
 1010 FORMAT(' Degeneracy for orbital ',I5,' of spin ',I3,' is ',I3)
 1020 FORMAT(I3,I3,I5,I5,F20.10,F20.10)
      END
