C compute relativistic corrections
      SUBROUTINE DRVREL(SCR,MAXCOR,NORBS)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION RELCOR(3)
      CHARACTER*6 LABELS(3)
      DIMENSION SCR(MAXCOR)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      DATA LABELS /'Darwin','p**4','Total'/
      I010=1+NORBS*NORBS
      CALL SEEKLB('DARWIN  ',IERR,0)
      IF(IERR.NE.0)RETURN
      CALL COMPPR(RELCOR(1), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('   P4   ',IERR,0)
      CALL COMPPR(RELCOR(2), SCR, SCR(I010),NORBS,.TRUE.)
      CALL SEEKLB('TOTAL   ',IERR,0)
      CALL COMPPR(RELCOR(3), SCR, SCR(I010),NORBS,.TRUE.)
      WRITE(LUOUT,100)
      WRITE(LUOUT,101)(LABELS(I),RELCOR(I),I=1,3)
100   FORMAT(T4,'Relativistic correction to the energy ')
101   FORMAT((T6,3(1X,A6, ' = ',F15.10)))
      RETURN
      END
