      SUBROUTINE SCLDEN(DENS,LDENS,SCR,LDIM,SCALE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION DENS(LDENS),SCR(LDIM,LDIM)
C
      CALL EXPND2(DENS,SCR,LDIM)
      DO 100 I=1,LDIM
        SCR(I,I)=SCALE*SCR(I,I)
  100 CONTINUE
      CALL SQUEZ2(SCR,DENS,LDIM)
      RETURN
      END
