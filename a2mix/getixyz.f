      SUBROUTINE GETIXYZ(IXYZ, IX, IY, IZ)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      IF (IXYZ .EQ. 1) THEN
         IX = 1
         IY = 0 
         IZ = 0
      ELSE IF (IXYZ .EQ. 2) THEN
         IX = 0
         IY = 1 
         IZ = 0
      ELSE 
         IX = 0
         IY = 0 
         IZ = 1
      ENDIF
C      
      RETURN
      END

