      Subroutine Energy_aa(T1aa,T1bb,T2aa,W2aa,Ndim_a,Ndim_b,Ncol,Nrow,
     +                     Listw,Irrep,Ispin,Tau,E,ES)

      Implicit Double Precision(A-H,O-Z)
      Logical tau

#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "info.com"
     
      Dimension T1aa(Ndim_a)
      Dimension T1bb(Ndim_b)
      Dimension T2aa(Nrow,Ncol)
      Dimension W2aa(Nrow,Ncol) 
      Dimension E(3),ES(3)
      Data One /1.0D0/

      Call Getlst(W2aa,1,Ncol,1,Irrep,Listw)

      If (Ispin .Eq. 1) Then
         ES(1) = ES(1) + Ddot(Nrow*Ncol,W2aa,1,T2aa,1)
      ELseif (Ispin .Eq. 2) then
         ES(2) = ES(2) + Ddot(Nrow*Ncol,W2aa,1,T2aa,1)
      Endif 

      If (Tau) then
         Call Ftau(T2aa,T1aa,T1bb,Nrow,Ncol,Pop(1,Ispin),
     +             Pop(1,Ispin),Vrt(1,Ispin),Vrt(1,Ispin),
     +             Irrep,Ispin,One) 
      Endif 

#ifdef _DEBUG_LVL0
      Print*, "After Tau"
      call checksum("Waa  :",W2aa,Ncol*Nrow)
      call checksum("T2aa :",T2aa,Ncol*Nrow)
#endif 

      If (Ispin .Eq. 1) Then
         E(1)  = E(1) + Ddot(Nrow*Ncol,W2aa,1,T2aa,1)
      ELseif (Ispin .Eq. 2) then
         E(2)  = E(2) + Ddot(Nrow*Ncol,W2aa,1,T2aa,1)
      Endif 
      
      Return
      End
      
