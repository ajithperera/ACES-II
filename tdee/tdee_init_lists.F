      Subroutine tdee_init_lists(work,Memleft,Iuhf,Irrepx,Source,
     +                           Target,Iside)

      Implicit Integer(A-Z)

      Dimension Work(Memleft)
      Dimension Target_lists(2),Source_lists(2)
      Logical Source, Target

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"
   
      Data Target_Lists /334,344/
      Data Source_Lists /314,324/

C Source and target singels

      pha_length = Irpdpd(Irrepx,9)
      phb_length = Irpdpd(Irrepx,10)
      hha_length = Irpdpd(Irrepx,21)
      hhb_length = Irpdpd(Irrepx,22)
      ppa_length = Irpdpd(Irrepx,19)
      ppb_length = Irpdpd(Irrepx,20)

      Max_s_length = Max(pha_length,phb_length,hha_length,hhb_length,
     +                   ppa_length,ppb_length)
     
      Call DZero(Work,Max_s_length)

      If (Source) Then
         Call Putlst(Work,1,1,1,1,390)
         Call Putlst(Work,1,1,1,1,392)
      Endif 

      If (Target) Then
         Call Putlst(Work,1,1,1,1,394)
         Call Putlst(Work,1,1,1,1,396)
         Call Putlst(Work,1,1,1,1,381)
         Call Putlst(Work,1,1,1,1,382)
      Endif 

      If (Iuhf .EQ. 1) Then
         If (Source) Then
            Call Putlst(Work,1,1,1,2,390)
            Call Putlst(Work,1,1,1,2,392)
         Endif

         If (Target) Then
            Call Putlst(Work,1,1,1,2,394)
            Call Putlst(Work,1,1,1,2,396)
            Call Putlst(Work,1,1,1,2,381)
            Call Putlst(Work,1,1,1,2,382)
         Endif 
      Endif 
      
C Source and target doubles lists

      If (Source) Then
         Do Ispin = 3, 3-2*Iuhf, -1
            Call Zerolist(Work,Memleft,(Source_lists(Iside)-1)+Ispin)
         Enddo 
      Elseif (Target) Then
         Do Ispin = 3, 3-2*Iuhf, -1
            Call Zerolist(Work,Memleft,(Target_lists(Iside)-1)+Ispin)
         Enddo 
      Endif

C  Other auxilary doubles lists

      If (Target) Then 

         If (Iuhf .EQ. 0) Then
            Call Zerolist(Work,Memleft,302)
         Else
            Call Zerolist(Work,Memleft,301)
            Call Zerolist(Work,Memleft,300)
            Call Zerolist(Work,Memleft,302) 
            Call Zerolist(Work,Memleft,305) 
            Call Zerolist(Work,Memleft,306) 
            Call Zerolist(Work,Memleft,308) 
         Endif 
         Call Zerolist(Work,Memleft,304)
         Call Zerolist(Work,Memleft,307)
         Call Zerolist(Work,Memleft,309)
         Call Zerolist(Work,Memleft,350)
         Call Zerolist(Work,Memleft,351)
         Call Zerolist(Work,Memleft,352)
         Call Zerolist(Work,Memleft,353)
         Call Zerolist(Work,Memleft,360)
         Call Zerolist(Work,Memleft,361)
         Call Zerolist(Work,Memleft,362)
         Call Zerolist(Work,Memleft,363)
      Endif 

      Return
      End

    
