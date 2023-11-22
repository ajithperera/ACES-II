      Subroutine Trans_vv_ab(T2_in,T2_out,Cvv,Nij,Irrepab,Iuhf,Type)

      Implicit Double Precision(A-H,O-Z)
      Integer Dissiz_trn,Dissiz_sqr
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
      Dimension T2_in(Nij),T2_out(Nij)

      Data Ione /1/
      Data Done,Dnull /1.0D0,0.0D0/ 

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

      Ioff = Ione
      Do Irrep = 1, Nirrep
         Irrep_b = Irrep
         Irrep_a = Dirprd(Irrep_b,Irrepab)
 
         Ncol = Vrt(Irrep_b,2)
         Nrow = Vrt(Irrep_a,1)

         Ioff_a = Ioff_vv(Irrep_a,1) 
         If (Iuhf .Eq. 0) Then
             Ioff_b = Ioff_vv(Irrep_b,1) 
         Else
             Ioff_b = Ioff_vv(Irrep_b,2) 
         Endif 

         Call Xgemm(A,B,Nrow,Ncol,Nrow,Done,Cvv(Ioff_a),Nrow,
     +              T2_in(Ioff),Nrow,Dnull,T2_out(Ioff),Nrow)

         Call Xgemm(C,D,Nrow,Ncol,Ncol,Done,T2_out(Ioff),Nrow,
     +              Cvv(Ioff_b),Ncol,Dnull,T2_in(Ioff),Nrow)

         Ioff = Ioff + Nrow*Ncol 
      Enddo

      Return
      End

 
