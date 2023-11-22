      Subroutine Check_source_lists(Work,Maxcor,Irrepx,Iuhf,Iside)
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

      Call Getlst(Work,1,1,1,1,390)
      Call checksum("list-390-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,2,390)
         Call checksum("list-390-BB      :",Work,Phb_length,s)
      Endif

      Elseif (Iside .Eq. 2) Then

      Call Getlst(Work,1,1,1,1,392)
      Call checksum("list-392-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,2,392)
         Call checksum("list-392-BB      :",Work,Phb_length,s)
      Endif 

      Endif 

      If (iside .eq. 1) Then 
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,314),isytyp(2,314))
      call getall(Work(1),length,irrepx,314)
      call checksum("List-314;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,315),isytyp(2,315))
      call getall(Work(1),length,irrepx,315)
      call checksum("List-315;Z(ab,ij):",Work(1),length,s)
      Endif
      Length=idsymsz(Irrepx,isytyp(1,316),isytyp(2,316))
      call getall(Work(1),length,irrepx,316)
      call checksum("List-316;Z(Ai,Bj):",Work(1),length,s)
      elseif (iside .eq. 2) then

      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,324),isytyp(2,324))
      call getall(Work(1),length,irrepx,324)
      call checksum("List-324;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,325),isytyp(2,325))
      call getall(Work(1),length,irrepx,325)
      call checksum("List-325;Z(ab,ij):",Work(1),length,s)
      Endif
      Length=idsymsz(Irrepx,isytyp(1,326),isytyp(2,326))
      call getall(Work(1),length,irrepx,326)
      call checksum("List-326;Z(Ai,Bj):",Work(1),length,s)
      Endif 

      Return 
      End
