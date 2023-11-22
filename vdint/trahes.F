      SUBROUTINE TRAHES(SHESS,CHESS,NCOOR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MATRIX
#include "mxcent.par"
      PARAMETER (HALF = 0.5D0)
      DIMENSION SHESS(NCOOR,NCOOR), CHESS(MXCOOR,MXCOOR),
     *          MATRIX(MXCOOR,MXCOOR)
      COMMON /TRANUC/ TRCTOS(MXCOOR,MXCOOR), TRSTOC(MXCOOR,MXCOOR)
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      DO 100 I = 1, NCOOR
         DO 200 J = 1, I - 1
            SHESS(J,I) = SHESS(I,J)
  200    CONTINUE
  100 CONTINUE
#ifdef _DEBUG_LVLM
      write(6,"(a)") " The symmetry Hessian"
      call output(SHESS,1,Ncoor,1,Ncoor,Ncoor,Ncoor,1,6)
#endif
      CALL XGEMM('N','N',NCOOR,NCOOR,NCOOR,ONE,
     *           TRSTOC,MXCOOR,SHESS,NCOOR,AZERO,
     *           MATRIX,MXCOOR)
      CALL XGEMM('N','T',NCOOR,NCOOR,NCOOR,ONE,
     *           MATRIX,MXCOOR,TRSTOC,MXCOOR,
     *           AZERO,CHESS,MXCOOR)
#ifdef _DEBUG_LVLM
      write(6,"(a)") " The no-symmetry Hessian"
      call output(CHESS,1,Ncoor,1,Ncoor,Mxcoor,Mxcoor,1,6)
#endif
      RETURN
      END
