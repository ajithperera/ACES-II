      Subroutine Pccd_print_all(Work,Nsize,Type,Irrepx)

      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Nsize) 
      Character*4 Type

#include "sym.com"
#include "sympop.com"
#include "syminf.com"

      Ioff = 1
      Do Irrep_r = 1, Nirrep
         Irrep_l = Dirprd(Irrep_r,Irrepx)
         If (Type .EQ. "ijkl") Then
            Nrow = Irpdpd(Irrep_l,21)
            Ncol = Irpdpd(Irrep_r,21)
            Call output(Work(Ioff),1,Nrow,1,Ncol,Nrow,Ncol,1)
	    Ioff = Ioff + Nrow*Ncol 
         Endif
         If (Type .EQ. "abcd") Then
            Nrow = Irpdpd(Irrep_l,19)
            Ncol = Irpdpd(Irrep_r,19)
            Call output(Work(Ioff),1,Nrow,1,Ncol,Nrow,Ncol,1)
	    Ioff = Ioff + Nrow*Ncol 
         Endif
      Enddo

      Return
      End 

      

