      Subroutine Pccd_reset_vv_ab(T2_in,Nab,Irrepab)

      Implicit Double Precision(A-H,O-Z)
      Integer Aoff

#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "symoff.com"
#include "sympop.com"
#include "info.com"

      Dimension T2_in(Nab)

      Data Ione,Dnull,Inull /1,0.0D0,0/

      Aoff = Ione 
      Do Irrep = 1, Nirrep

         Irrep_b = Irrep
         Irrep_a = Dirprd(Irrep_b,Irrepab)
 
         Ncol = Vrt(Irrep_b,1)
         Nrow = Vrt(Irrep_a,1)
         If (Irrep_a .Eq. Irrep_b) Then
            Call Pccd_set2_zero(T2_in(Aoff),Nrow,Ncol)
            Aoff = Aoff + Nrow*Ncol
         Else
            Call Dzero(T2_in(Aoff),Nrow*Ncol)
            Aoff = Aoff + Nrow*Ncol
         Endif 

      Enddo

      Return
      End

 
