      SUBROUTINE TRAHES(SHESS,CHESS,MATRIX,TRSTOC,NCOOR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MATRIX
C
      DIMENSION SHESS(NCOOR,NCOOR), CHESS(NCOOR,NCOOR),
     &          MATRIX(NCOOR,NCOOR),TRSTOC(NCOOR,NCOOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DATA AZERO,ONE /0.D0,1.D0/
C
C READ IN THE TRANSFORMATION MATRIX
C
      CALL GETREC(20,'JOBARC','SYMCOORD',IINTFP*NCOOR*NCOOR,
     &            TRSTOC)
C
C TRANSFORM RIGHT AND LEFT HAND SIDE COORDINATES
C
C
      CALL XGEMM('N','N',NCOOR,NCOOR,NCOOR,ONE,TRSTOC,NCOOR,
     &           SHESS,NCOOR,AZERO,MATRIX,NCOOR)
      CALL XGEMM('N','T',NCOOR,NCOOR,NCOOR,ONE,MATRIX,NCOOR,
     &           TRSTOC,NCOOR,AZERO,CHESS,NCOOR)
C
      RETURN
      END
