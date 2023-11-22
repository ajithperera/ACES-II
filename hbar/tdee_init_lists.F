      Subroutine tdee_init_lists(work,Memleft,Iuhf)

      Implicit Integer(A-Z)

      Dimension Work(Memleft)
      Dimension Target_lists(2),Source_lists(2)
      Logical Source, Target

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"
   
      Irrepx = 1

      pha_length = Irpdpd(Irrepx,9)
      phb_length = Irpdpd(Irrepx,10)
      hha_length = Irpdpd(Irrepx,21)
      hhb_length = Irpdpd(Irrepx,22)
      ppa_length = Irpdpd(Irrepx,19)
      ppb_length = Irpdpd(Irrepx,20)

      Max_s_length = Max(pha_length,phb_length,hha_length,hhb_length,
     +                   ppa_length,ppb_length)
     
      Call DZero(Work,Max_s_length)

      Call Putlst(Work,1,1,1,1,90)
      Call Putlst(Work,1,1,1,2,90)

C  Other auxilary doubles lists

      If (Iuhf .EQ. 0) Then
        Call Zerolist(Work,Memleft,34)
        Call Zerolist(Work,Memleft,37)
        Call Zerolist(Work,Memleft,39)
        Call Zerolist(Work,Memleft,42)
        Call Zerolist(Work,Memleft,43)
        Call Zerolist(Work,Memleft,44)
        Call Zerolist(Work,Memleft,46)
      Else
        Call Zerolist(Work,Memleft,34)
        Call Zerolist(Work,Memleft,35)
        Call Zerolist(Work,Memleft,36)
        Call Zerolist(Work,Memleft,37)
        Call Zerolist(Work,Memleft,38)
        Call Zerolist(Work,Memleft,39)
        Call Zerolist(Work,Memleft,40)
        Call Zerolist(Work,Memleft,41)
        Call Zerolist(Work,Memleft,42)
        Call Zerolist(Work,Memleft,43)
        Call Zerolist(Work,Memleft,44)
        Call Zerolist(Work,Memleft,45)
        Call Zerolist(Work,Memleft,46)
      Endif 

      Return
      End

    
