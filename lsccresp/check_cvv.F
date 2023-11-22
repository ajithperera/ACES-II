      Subroutine check_cvv(Cvv,Iuhf)

      Implicit Double Precision (A-H,O-Z)

#include "syminf.com"
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
#include "info.com"
      
      Dimension Cvv(Nfea(1)+Iuhf*Nfea(2))
     
      Do Ispin = 1, Iuhf+1 
         Ibegin = 1 + Nfea(1)*(Ispin-1)
         Ioff = 0
         Write(6,"(a)") " The CVV matrix"
      Do Irrep = 1, Nirrep
         Ndim = Vrt(Irrep,Ispin) 
         Call output(Cvv(Ioff+Ibegin),1,Ndim,1,Ndim,Ndim,Ndim,1)
         Ioff = Ioff + Ndim*Ndim
      Enddo 
      Enddo

      Return
      End 
