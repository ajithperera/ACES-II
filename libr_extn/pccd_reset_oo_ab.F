      Subroutine Pccd_reset_oo_ab(T2_in,Nij,Irrepab)

      Implicit Double Precision(A-H,O-Z)

#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "symoff.com"
#include "sympop.com"
#include "info.com"

      Dimension T2_in(Nij)

      Data Ione,Dnull,Inull /1,0.0D0,0/

      Ioff = Ione
      Do Irrep = 1, Nirrep
         Irrep_b = Irrep
         Irrep_a = Dirprd(Irrep_b,Irrepab)
         Ncol = Pop(Irrep_b,1)
         Nrow = Pop(Irrep_a,1)
         If (Irrep_a .Eq. Irrep_b) Then 
             Call Pccd_set2_zero(T2_in(Ioff),Nrow,Ncol)
             Ioff = Ioff + Nrow*Ncol
         Else
             Call Dzero(T2_in(Ioff),Nrow*Ncol)
             Ioff = Ioff + Nrow*Ncol
         Endif 
      Enddo

      Return
      End

 
