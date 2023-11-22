      Subroutine Zero_490_lists(Work,Length,Irrepx,Iuhf)
 
      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Length)

      Integer pha_length,phb_length

#include "sym.com"
#include "sympop.com"
#include "syminf.com"
  
      Pha_length = Irpdpd(Irrepx,9)
      Phb_length = Irpdpd(Irrepx,10)

      Max_length = Max(Pha_length,Phb_length)

      Call Dzero(Work,Max_length)
      Call Putlst(Work,1,1,1,3,490)

      If (Iuhf .Ne. 0) Then
         Call Putlst(Work,1,1,1,4,490)
      Endif

      Return
      End
