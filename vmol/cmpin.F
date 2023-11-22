      SUBROUTINE CMPIN(A,AC,TOL,N,IH,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "baslims.par"
      parameter (nh4=4*nht-3, nh2=nht+nht+1)
      parameter (mxp2=maxprim*maxprim)
      parameter (mxp4=mxp2*mxp2)
      DIMENSION A(2),B(2),AC(2),BC(2),T(2)
      COMMON /ICOMP/ IBE(mxp4),K,NN
      common /icomp2/IBEM
      NN=N
      DO 1 I=1,N
    1   T(I) = ABS(A(I))
      CALL WHENFGT(N,T,1,TOL,IBE,IBEM)
      CALL GATHER(IBEM,AC,A,IBE)
      IH=IBEM
      RETURN
      ENTRY CMPR1(A,AC)
      CALL GATHER(IBEM,AC,A,IBE)
      RETURN
      ENTRY CMPR2(A,B,AC,BC)

      CALL GATHER(IBEM,AC,A,IBE)
      CALL GATHER(IBEM,BC,B,IBE)
      RETURN
      ENTRY EXPND(A,AC)
      CALL SCATTER(IBEM,A,IBE,AC)
      RETURN
      ENTRY CMPRI(A,AC,M)
      DO 44 MM=0,M-1
   44   CALL GATHER(IBEM,AC(MM*ibem+1),A(MM*NN+1),IBE)
      RETURN
      END
