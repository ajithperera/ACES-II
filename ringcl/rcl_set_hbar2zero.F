      Subroutine Rcl_set_hbar2zero(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"

      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)

C This routine selectively set Hbar elements to zero.

C Hbar(MN,IJ)
      IF (UHF) Then
         Call Zerolist(Work,Maxcor,51)
         Call Zerolist(Work,Maxcor,52)
         Call Zerolist(Work,Maxcor,53)
      Else
         Call Zerolist(Work,Maxcor,53)
      Endif

C These are needed from doubles correction.
C Hbar(AI,BC)
CSSS      IF (UHF) Then
CSSS         Call Zerolist(Work,Maxcor,27)
CSSS         Call Zerolist(Work,Maxcor,28)
CSSS         Call Zerolist(Work,Maxcor,29)
CSSS         Call Zerolist(Work,Maxcor,30)
CSSS      Else
CSSS         Call Zerolist(Work,Maxcor,27)
CSSS         Call Zerolist(Work,Maxcor,30)
CSSS      Endif 

C Hbar(IJ,KA)
CSSS      IF (UHF) Then
CSSS         Call Zerolist(Work,Maxcor,7)
CSSS         Call Zerolist(Work,Maxcor,8)
CSSS         Call Zerolist(Work,Maxcor,9)
CSSS         Call Zerolist(Work,Maxcor,10)
CSSS      Else
CSSS         Call Zerolist(Work,Maxcor,7)
CSSS         Call Zerolist(Work,Maxcor,10)
CSSS      Endif

C Hbar(AB,CI)
      IF (UHF) Then
         Call Zerolist(Work,Maxcor,127)
         Call Zerolist(Work,Maxcor,128)
         Call Zerolist(Work,Maxcor,129)
         Call Zerolist(Work,Maxcor,130)
      Else
         Call Zerolist(Work,Maxcor,130)
      Endif

C Hbar(IA,JK)
      IF (UHF) Then
         Call Zerolist(Work,Maxcor,107)
         Call Zerolist(Work,Maxcor,108)
         Call Zerolist(Work,Maxcor,109)
         Call Zerolist(Work,Maxcor,110)
      Else
CSSS         Call Zerolist(Work,Maxcor,107)
         Call Zerolist(Work,Maxcor,110)
      Endif

      Return
      End

