      subroutine Trans_oo_ab(T2_in,T2_out,Coo,Nij,Irrepij,Iuhf,Type)

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

      Dimension Coo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension T2_in(Nij),T2_out(Nij)

      Data Done,Dnull /1.0D0,0.0D0/
      Data Ione /1/

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

      Ioff   = Ione
      Do Irrep = 1, Nirrep
         Irrep_j = Irrep
         Irrep_i = Dirprd(Irrep_j,Irrepij)
 
         Nrow = Pop(Irrep_i,1)
         Ncol = Pop(Irrep_j,2)

         Ioff_a = Ioff_oo(Irrep_i,1) 
         If (Iuhf .Eq. 0) Then
            Ioff_b = Ioff_oo(Irrep_j,1)
         Else
            Ioff_b = Ioff_oo(Irrep_j,2) 
         Endif 

         Call Xgemm(A,B,Nrow,Ncol,Nrow,Done,Coo(Ioff_a),Nrow,
     +              T2_in(Ioff),Nrow,Dnull,T2_out(Ioff),Nrow)

         Call Xgemm(C,D,Nrow,Ncol,Ncol,Done,T2_out(Ioff),Nrow,
     +              Coo(Ioff_b),Ncol,Dnull,T2_in(Ioff),Nrow)

         Ioff  = Ioff + Nrow*Ncol
      Enddo

      Return
      End

 
