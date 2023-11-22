      Subroutine Rcc_dbg(Work,Maxcor,Iuhf,Ncycle)
 
      Implicit Double Precision(A-H,O-Z)
      Dimension Work(Maxcor)

#include "flags.com"
#include "machsp.com"
#include "sympop.com"
#include "syminf.com"
  
      If (Ncycle .EQ. 1) Then
      Write(6,"(a,a)") " Ringcl debugging is on and F(ab) and F(ij)",
     +                 " is added."
      Write(6,"(a)") " These terms are not naturally in RCC or DRCC "
      Write(6,"(a)") " ------????????????-----------"
      Endif 

      Call Rcc_fab(Work,Maxcor,Iuhf,1.0D0)
      Call Rcc_fij(Work,Maxcor,Iuhf,1.0D0)

      Call Rcc_fab_int2(Work,Maxcor,Iuhf)
      Call Rcc_fij_int2(Work,Maxcor,Iuhf)

      Return
      End 
