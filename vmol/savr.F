      SUBROUTINE SAVR (AREA,NAQRS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Parameters
C-----------------------------------------------------------------------
#include "baslims.par"
      parameter (nh4=4*nht-3, nh2=nht+nht+1)
      parameter (mxp2=maxprim*maxprim)
C-----------------------------------------------------------------------
      DIMENSION AREA(NAQRS,4)
      COMMON /ONE/ 
     1 ALP(mxp2),BET(mxp2),SS1(mxp2),CSS1(mxp2),PV(mxp2,8,3),CAB(mxp2),
     2 GAM(mxp2),DEL(mxp2),SS2(mxp2),CSS2(mxp2),QV(mxp2,8,3),CCD(mxp2)
      COMMON /CON/ RAX,RAY,RAZ,RPRX,RPRY,RPRZ,DISTPQ,DISTRS,TOLR,TP52,
     &             ZJUNK(MXP2,16)
      COMMON /CONI/NUCT(4),NRCT(4),JSTT(4),
     1 IFD1X,IFD2X,NUCAQ,NUCRS,NXXXX(2),KQ,KR,KS,IVA
      CALL MVSUM  (AREA(1,3),SS1,NUCAQ,SS2,NUCRS)
      CALL MVPROD (AREA(1,4),SS1,NUCAQ,SS2,NUCRS)
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 2 I = 1,NAQRS
C
      AREA(I,3) = 1./AREA(I,3)
    2 AREA(I,2) = AREA(I,4) * AREA(I,3)
      CALL MVPROD (AREA(1,4),CSS1,NUCAQ,CSS2,NUCRS)
CDIR$ IVDEP
*VOCL LOOP,NOVREC
      DO 3 I=1,NAQRS
      AREA(I,1)= AREA(I,4)*SQRT(AREA(I,3))*TP52
    3 CONTINUE
      RETURN
      END
