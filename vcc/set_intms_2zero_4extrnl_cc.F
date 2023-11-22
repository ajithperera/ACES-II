      Subroutine Set_intms_2zero_4extrnl_cc(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"

      UHF    = (IUhf .NE. 0)
      IRREPX = 1 

C First Hbar(A,E), Hbar(M,I) and Hbar(M,E) 

      HHA_LENGTH = IRPDPD(IRREPX,21)
      HHB_LENGTH = IRPDPD(IRREPX,22)

      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)

      PPA_LENGTH = IRPDPD(IRREPX,19)
      PPB_LENGTH = IRPDPD(IRREPX,20)

      Call Zerlst(Work,HHA_LENGTH,1,1,1,91)
      Call Zerlst(Work,PPA_LENGTH,1,1,1,92)
      Call Zerlst(Work,HPA_LENGTH,1,1,1,93)

      If (Uhf) Then
         Call Zerlst(Work,HHB_LENGTH,1,1,2,91)
         Call Zerlst(Work,PPB_LENGTH,1,1,2,92)
         Call Zerlst(Work,HPB_LENGTH,1,1,2,93)
       Endif 

      Call Zerolist(Work,Maxcor,53)
      Call Zerolist(Work,Maxcor,54)
      Call Zerolist(Work,Maxcor,56)
      Call Zerolist(Work,Maxcor,58)
      Call Zerolist(Work,Maxcor,61)
      Call Zerolist(Work,Maxcor,63)

      If (Uhf) Then
         Call Zerolist(Work,Maxcor,51)
         Call Zerolist(Work,Maxcor,52)
         Call Zerolist(Work,Maxcor,62)
         Call Zerolist(Work,Maxcor,55)
         Call Zerolist(Work,Maxcor,57)
         Call Zerolist(Work,Maxcor,59)
      Endif 

      Return
      End


