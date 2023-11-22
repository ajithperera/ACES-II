      FUNCTION DERSTP(GDDIA,HESDIA,RNU,NCORD,DUMMY)
C
C     Purpose:
C
C        Calculate derivative of step length function with
C        respect to RNU
C
C     DERSTP = STEP*STEP' / //STEP//
C
C        where
C
C     STEP = -GDDIA / (HESDIA+RNU)
C     STEP'=  GDDIA / ((HESDIA+RNU)**2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GDDIA(2),HESDIA(2)
      PARAMETER ( D0=0.0D0 )
      DERSTP=D0
      DEL   =D0
      DO 100 K = 1,NCORD
         DEL    = DEL    + (GDDIA(K)/(HESDIA(K)+RNU))**2
         DERSTP = DERSTP - (GDDIA(K)**2) / ((HESDIA(K)+RNU)**3)
 100  CONTINUE
      DEL    = SQRT(DEL)
      DERSTP = DERSTP/DEL
      RETURN
      END
