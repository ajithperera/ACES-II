      Subroutine Pccd_transp_pphh(Pphh,Hhpp,Nsize)

      Implicit Double Precision(A-H,O-Z)
      Dimension Pphh(Nsize)
      Dimension Hhpp(Nsize)

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Data One,Onem,Dnull,Ione /1.0D0,-1.0D0,0.0D0,1/

      Irrepx = 1

C H(ab,ij) -> H(ij,ab)

      Ioff = Ione 
      Do Irrep_ij = 1, Nirrep
         Irrep_ab = Dirprd(Irrep_ij,Irrepx)

         Nrow_ab = Irpdpd(Irrep_ab,15)
         Ncol_ij = Irpdpd(Irrep_ij,14)
         Call Transp(Pphh(Ioff),Hhpp(Ioff),Ncol_ij,Nrow_ab)
         Ioff = Ioff + Nrow_ab*Ncol_ij

      Enddo 

      Return
      End
