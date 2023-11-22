#include "flags.h"
      Subroutine Rcc_expand(Win,Wout,Ncol_squar_ij,Nrow_squar_ab,
     +                      ncol_tring_ij,Nrow_tring_ab,
     +                      Irrep,Ispin,Work,Left)

      Implicit Double Precision (A-H, O-Z)
      Dimension Work(Left),Win(ncol_tring_ij*Nrow_tring_ab)
      Dimension Wout(Ncol_squar_ij*Nrow_squar_ab)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "parcc.com"
#include "flags.com"
#include "info.com"

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") "--------Entering rcc_expand---------"
      Write(6,*)
#endif 
      Irrepl = Irrep
      Irrepr = Irrep

      I000 = 1
      I010 = I000 + Nrow_squar_ab*Ncol_squar_ij
      I020 = I010 + Nrow_squar_ab*Ncol_squar_ij
      Iend = I020
      IF (Iend .Gt. Left) Call Insmem("rcc_expand",Iend,Left)

C D(A<=B,I<=J) -->D(AB,I<=J)

      Call Symexp6(Irrepl,Vrt(1,ispin),Vrt(1,ispin),
     +             Nrow_squar_ab,Nrow_tring_ab,
     +             ncol_tring_ij,work(I000),Win,
     +             work(I010))
      Call Transp(work(I000),work(I010),ncol_tring_ij,
     +            nrow_squar_ab)
C D(IJ,A<=B) -->D(iJ,AB)

      Call Symexp6(Irrepr,Pop(1,ispin),Pop(1,ispin),
     +             Ncol_squar_ij,Ncol_tring_ij,
     +             nrow_squar_ab,work(I000),work(I010),
     +             Wout)
C D(IJ,AB) -->D(AB,IJ)

       Call Transp(work(I000),Wout,Nrow_squar_ab,
     +                  Ncol_squar_ij)
#ifdef _DEBUG_LVLM
      Write(6,*) "D<AB|IJ>"
      call output(Wout,1,Nrow_squar_ab,1,Ncol_squar_ij,
     +            Nrow_squar_ab,Ncol_squar_ij,1)
#endif 
      Return
      End

