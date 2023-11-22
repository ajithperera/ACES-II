      Subroutine Check_target_lists(Work,Maxcor,Irrepx,Iuhf,Iside)
      Implicit Integer(A-Z)

      Double Precision Work(Maxcor)
      Logical UHF

#include "machsp.com"
#include "sympop.com"

      pha_length = Irpdpd(Irrepx,9)
      phb_length = Irpdpd(Irrepx,10)
      hha_length = Irpdpd(Irrepx,21)
      hhb_length = Irpdpd(Irrepx,22)
      ppa_length = Irpdpd(Irrepx,19)
      ppb_length = Irpdpd(Irrepx,20)

      Write(6,*) 

      If (iside .Eq. 1) then 

      Call Getlst(Work,1,1,1,1,394)
      Call checksum("list-394-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,2,394)
         Call checksum("list-394-BB      :",Work,Phb_length,s)
      Endif

      Elseif (Iside .Eq. 2) Then

      Call Getlst(Work,1,1,1,1,396)
      Call checksum("list-396-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,2,396)
         Call checksum("list-396-BB      :",Work,Phb_length,s)
      Endif 

      Endif 

      If (iside .eq. 1) Then 
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,334),isytyp(2,334))
      call getall(Work(1),length,irrepx,334)
      call checksum("List-334;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,335),isytyp(2,335))
      call getall(Work(1),length,irrepx,335)
      call checksum("List-335;Z(ab,ij):",Work(1),length,s)
      Endif
      Length=idsymsz(Irrepx,isytyp(1,336),isytyp(2,336))
      call getall(Work(1),length,irrepx,336)
      call checksum("List-336;Z(Ai,Bj):",Work(1),length,s)
      elseif (iside .eq. 2) then

      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,344),isytyp(2,344))
      call getall(Work(1),length,irrepx,344)
      call checksum("List-344;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,335),isytyp(2,335))
      call getall(Work(1),length,irrepx,345)
      call checksum("List-345;Z(ab,ij):",Work(1),length,s)
      Endif
      Length=idsymsz(Irrepx,isytyp(1,346),isytyp(2,346))
      call getall(Work(1),length,irrepx,346)
      call checksum("List-346;Z(Ai,Bj):",Work(1),length,s)
      Endif 

      Return 
      End
