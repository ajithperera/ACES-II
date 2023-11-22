      SUBROUTINE HARMTR(T,MAXL,E)
      IMPLICIT NONE
C-----------------------------------------------------------------------
      DOUBLE PRECISION T,E
      INTEGER MAXL
C-----------------------------------------------------------------------
      DOUBLE PRECISION COST,SINT,COSP,SINP
C-----------------------------------------------------------------------
      DIMENSION T(9,9),E(3)
C
      write(6,*) ' e vector ',e
C
      COST = E(3)
      IF( (1.0D+00-COST**2) .GT. 1.0D-10)THEN
       SINT = DSQRT( 1.0D+00-COST**2 )
      ELSE
       SINT = 0.0D+00
      ENDIF
C
      IF(SINT .GT. 1.0D-06)THEN
       COSP = E(1)/SINT
       SINP = E(2)/SINT
      ELSE
       COSP = 1.0D+00
       SINP = 0.0D+00
      ENDIF
C
      write(6,*) ' cosp, sinp, cost, sint '
      write(6,*) cosp,sinp,cost,sint
C
      CALL ZERO(T,9*9)
      T(1,1) = 1.0D+00
      T(2,2) = 1.0D+00
      T(3,3) = 1.0D+00
      T(4,4) = 1.0D+00
c      return
C
      IF(MAXL .EQ. 0)THEN
       write(6,*) ' t matrix '
       call output(t,1,9,1,9,9,9,1)
       RETURN
      ELSEIF(MAXL.EQ.1)THEN
       T(2,2) =  COST*COSP
       T(2,3) = -SINP
       T(2,4) =  SINT*COSP
       T(3,2) =  COST*SINP
       T(3,3) =  COSP
       T(3,4) =  SINT*SINP
       T(4,2) = -SINT
       T(4,4) =  COST
       write(6,*) ' t matrix '
       call output(t,1,9,1,9,9,9,1)
      ENDIF
C
      RETURN
      END
