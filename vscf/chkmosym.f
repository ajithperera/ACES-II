      SUBROUTINE CHKMOSYM(IRRSYM,NBAS,IUHF,MOSYMOK)
C
C  This routine reads records EVCSYMAC and EVCSYMBC to see if the AO
C  basis MO coefficients transform as irreducible representations of
C  the computational point group.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NBAS,IUHF
      CHARACTER*8 IRRSYM(NBAS)
      CHARACTER*4 PTGRP
      LOGICAL MOSYMOK
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C     Check if MOs transform as irreducible representations of the
C     computational point group.
C
      CALL GETCREC(20,'JOBARC','COMPPTGP',4,PTGRP)
C
      MOSYMOK = .TRUE.
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
      IF( IRRSYM(IMO)(1:4) .EQ. '    ' .OR.
     &    IRRSYM(IMO)(1:4) .EQ. '   u' .OR.
     &    IRRSYM(IMO)(1:4) .EQ. '   g'       )THEN
       WRITE(6,1000) ISPIN,IMO
       MOSYMOK = .FALSE.
      ENDIF
C
C     Special Cs group code.
C
      IF(PTGRP.EQ.'C s ' .OR. PTGRP.EQ.'C1h ')THEN
       IF(IRRSYM(IMO)(1:4) .EQ. ' A  ')THEN
        WRITE(6,1000) ISPIN,IMO
        WRITE(6,'(A4)') IRRSYM(IMO)(1:4)
        MOSYMOK = .FALSE.
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
        MOSYMOK = .FALSE.
       ENDIF
      ENDIF
C
   10 CONTINUE
   20 CONTINUE
      RETURN
 1000 FORMAT(' @CHKMOSYM-I, Spin ',I3,' Orbital ',I5,
     &       ' is not symmetry-adapted. ')
      END
