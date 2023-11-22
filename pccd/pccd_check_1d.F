      Subroutine Pccd_check_1d(Label,Hpq,N,Type,S)

      Implicit Double Precision(A-H,O-Z)
      Character*9 Label 
      Character*2 Type
      Character*1 S
      Dimension Hpq(N*N)

#include "syminf.com"   
#include "sym.com"   

      Data Dnull,Ione,Inull /0.0D0,1,0/

      E   = Dnull
      Ioff= Ione
      do Irrep = 1, Nirrep
         If (Type .Eq. "OO") M = Pop(Irrep,1)
         If (Type .Eq. "VV") M = Vrt(Irrep,1)
         Call Pccd_sum(Hpq(Ioff),M,E,S)
         Ioff = Ioff + M*M
      enddo

      write(6,"(a,2(2x,F15.10))") Label, e,e**2 
 
      return
      end 
