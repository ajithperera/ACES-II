      SUBROUTINE ERNGAB(T,Z,NSIZE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B
      DIMENSION T(NSIZE),Z(NSIZE)
      COMMON /FLAGS/ IFLAGS(100)
      CALL GETALL(T,NSIZE,1,46)
      X=SDOT(NSIZE,Z,1,T,1)
      IF(IFLAGS(1).GE.10)THEN
       WRITE(6,310)X
      ENDIF
310   FORMAT(T3,'W(mbej) AB contribution =',F14.10,' a.u.')
      RETURN
      END
