      Subroutine Check_lists(Work,Maxcor,Irrepx,Iuhf,Iside)
      Implicit Integer(A-Z)

      Double Precision Work(Maxcor)
      Logical UHF

#include "machsp.com"
#include "sympop.com"

#ifdef _DEBUG_LVL0 
      Write(6,*)
      Length=idsymsz(Irrepx,isytyp(1,304),isytyp(2,304))
      call getall(Work(1),length,irrepx,304)
      call checksum("List-304",Work(1),length,s)
      If (iuhf .gt.0) Then
      Length=idsymsz(Irrepx,isytyp(1,305),isytyp(2,305))
      call getall(Work(1),length,irrepx,305)
      call checksum("List-305",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,306),isytyp(2,306))
      call getall(Work(1),length,irrepx,306)
      call checksum("List-306",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,308),isytyp(2,308))
      call getall(Work(1),length,irrepx,308)
      call checksum("List-308",Work(1),length,s)
      endif 
      Length=idsymsz(Irrepx,isytyp(1,307),isytyp(2,307))
      call getall(Work(1),length,irrepx,307)
      call checksum("List-307",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,309),isytyp(2,309))
      call getall(Work(1),length,irrepx,309)
      call checksum("List-309",Work(1),length,s)

      Length=idsymsz(Irrepx,isytyp(1,350),isytyp(2,350))
      call getall(Work(1),length,irrepx,350)
      call checksum("List-350",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,351),isytyp(2,351))
      call getall(Work(1),length,irrepx,351)
      call checksum("List-351",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,352),isytyp(2,352))
      call getall(Work(1),length,irrepx,352)
      call checksum("List-352",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,353),isytyp(2,353))
      call getall(Work(1),length,irrepx,353)
      call checksum("List-353",Work(1),length,s)

      If (iside .eq. 2) Then
      Length=idsymsz(Irrepx,isytyp(1,360),isytyp(2,360))
      call getall(Work(1),length,irrepx,360)
      call checksum("List-360",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,361),isytyp(2,361))
      call getall(Work(1),length,irrepx,361)
      call checksum("List-361",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,362),isytyp(2,362))
      call getall(Work(1),length,irrepx,362)
      call checksum("List-362",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,363),isytyp(2,363))
      call getall(Work(1),length,irrepx,363)
      Endif 
#endif 
#ifdef _DEBUG_LVLM
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,334),isytyp(2,334))
      call getall(Work(1),length,irrepx,334)
      call checksum("List-334",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,335),isytyp(2,335))
      call getall(Work(1),length,irrepx,335)
      call checksum("List-335",Work(1),length,s)
      Endif
      Length=idsymsz(Irrepx,isytyp(1,336),isytyp(2,336))
      call getall(Work(1),length,irrepx,336)
      call checksum("List-336",Work(1),length,s)
      If (iside .eq. 2) Then 
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,344),isytyp(2,344))
      call getall(Work(1),length,irrepx,344)
      call checksum("List-344",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,335),isytyp(2,335))
      call getall(Work(1),length,irrepx,345)
      call checksum("List-345",Work(1),length,s)
      Endif
      Length=idsymsz(Irrepx,isytyp(1,346),isytyp(2,346))
      call getall(Work(1),length,irrepx,346)
      call checksum("List-346",Work(1),length,s)
      Endif 
#endif 
      Return 
      End
