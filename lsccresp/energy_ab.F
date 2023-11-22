      Subroutine Energy_ab(T1aa,T1bb,T2ab,W2ab,Ndim_a,Ndim_b,Nrow,Ncol,
     +                     Listw,Irrep,Tau,E,ES)

      Implicit Double Precision(A-H,O-Z)
      Logical tau

#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
     
      Dimension T1aa(Ndim_a)
      Dimension T1bb(Ndim_b)
      Dimension T2ab(Nrow,Ncol)
      Dimension W2ab(Nrow,Ncol) 
      Dimension E(3),ES(3)
      Data One /1.0D0/

      Call Getlst(W2ab,1,Ncol,1,Irrep,Listw)

      ES(3) = ES(3) + Ddot(Nrow*Ncol,W2ab,1,T2ab,1)

      If (Tau) then
         Call Ftau(T2ab,T1aa,T1bb,Nrow,Ncol,Pop(1,1),Pop(1,2),
     +             Vrt(1,1),Vrt(1,2),Irrep,3,One)
      Endif 

#ifdef _DEBUG_LVL0
      call checksum("Wab  :",W2ab,Ncol*Nrow)
      call checksum("T2ab :",T2ab,Ncol*Nrow)
#endif
      E(3)  = E(3) + Ddot(Nrow*Ncol,W2ab,1,T2ab,1)

      Return
      End
      
