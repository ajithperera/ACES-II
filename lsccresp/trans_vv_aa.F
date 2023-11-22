      Subroutine Trans_vv_aa(T2_in,T2_out,Cvv,Nab,Irrepab,Iuhf,Ispin,
     +                       Type)

      Implicit Double Precision(A-H,O-Z)
      Character*1 A,B,C,D
      Character*3 Type 

#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "symoff.com"
#include "sympop.com"
#include "info.com"

      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension T2_in(Nab),T2_out(Nab)

      Data Ione /1/
      Data Done,Dnull /1.0D0,0.0D0/ 

#ifdef _DEBUG_LVLM
      call checksum("T2vv-I:",T2_in,Nab)
#endif

      If (Type .Eq. "C2N") Then
          A = "T"
          B = "N"
          C = "N"
          D = "N"
      elseif (Type .Eq. "N2C") Then
          A = "N"
          B = "N"
          C = "N"
          D = "T"
      Endif 

      Koff = Ione
      Do Irrep = 1, Nirrep
         Irrep_b = Irrep
         Irrep_a = Dirprd(Irrep_b,Irrepab)
 
         Ncol = Vrt(Irrep_b,Ispin)
         Nrow = Vrt(Irrep_a,Ispin)

         Ioff = Ioff_vv(Irrep_a,Ispin) 
         Joff = Ioff_vv(Irrep_b,Ispin) 

         Call Xgemm(A,B,Nrow,Ncol,Nrow,Done,Cvv(Ioff),Nrow,
     +              T2_in(Koff),Nrow,Dnull,T2_out(Koff),Nrow)

         Call Xgemm(C,D,Nrow,Ncol,Ncol,Done,T2_out(Koff),Nrow,
     +              Cvv(Joff),Ncol,Dnull,T2_in(Koff),Nrow)

         Koff = Koff + Nrow*Ncol
      Enddo

#ifdef _DEBUG_LVLM
      call checksum("T2vv-O:",T2_in,Nab)
#endif
      Return
      End

 
