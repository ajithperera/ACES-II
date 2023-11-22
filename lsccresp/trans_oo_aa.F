      Subroutine Trans_oo_aa(T2_in,T2_out,Coo,Nij,Irrepij,Iuhf,Ispin,
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

      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension T2_in(Nij),T2_out(Nij)

      Data Ione /1/
      Data Done,Dnull /1.0D0,0.0D0/

#ifdef _DEBUG_LVLM
      call checksum("T2-in:",T2_in,Nij)
#endif
C Canonical to NO transformation.

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
         Irrep_j = Irrep
         Irrep_i = Dirprd(Irrep_j,Irrepij)
 
         Ncol = Pop(Irrep_j,Ispin)
         Nrow = Pop(Irrep_i,Ispin)

         Ioff = Ioff_oo(Irrep_i,Ispin) 
         Joff = Ioff_oo(Irrep_j,Ispin) 

         Call Xgemm(A,B,Nrow,Ncol,Nrow,Done,Coo(Ioff),Nrow,
     +              T2_in(Koff),Nrow,Dnull,T2_out(Koff),Nrow)

         Call Xgemm(C,D,Nrow,Ncol,Ncol,Done,T2_out(Koff),Nrow,
     +              Coo(Joff),Ncol,Dnull,T2_in(Koff),Nrow)

         Koff = Koff + Nrow*Ncol
      Enddo

#ifdef _DEBUG_LVLM
      call checksum("T2out:",T2_in,Nij)
#endif
      Return
      End

 
