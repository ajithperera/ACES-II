      FUNCTION FUZZY(N,F0)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DO I=1,N
C         F1=1.5D+00*F0-0.5D+00*F0**3.D+00
C
C Ajith 05/29/96 DEC Alpha do not like this
C
         F1=1.5D+00*F0-0.5D+00*F0**3
         F0=F1
      END DO
      FUZZY=F0
C
      RETURN
      END
