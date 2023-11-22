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

      Call Getlst(Work,1,1,1,3,490)
      Call checksum("list-490-AA      :",Work,Pha_length)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,4,490)
         Call checksum("list-490-BB      :",Work,Phb_length)
      Endif

      Elseif (Iside .Eq. 2) Then

      Call Getlst(Work,1,1,1,3,490)
      Call checksum("list-490-AA      :",Work,Pha_length)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,4,490)
         Call checksum("list-490-BB      :",Work,Phb_length)
      Endif 

      Endif 

      If (iside .eq. 1) Then 
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,461),isytyp(2,461))
      call getall(Work(1),length,irrepx,461)
      call checksum("List-461;Z(AB,IJ):",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,462),isytyp(2,462))
      call getall(Work(1),length,irrepx,462)
      call checksum("List-462;Z(ab,ij):",Work(1),length)
      Endif
      Length=idsymsz(Irrepx,isytyp(1,463),isytyp(2,463))
      call getall(Work(1),length,irrepx,463)
      call checksum("List-463;Z(Ai,Bj):",Work(1),length)

C check the temp. ring lists.
#ifdef _NOSKIP
      If (iuhf .ne.0) then
      Length=idsymsz(Irrepx,isytyp(1,434),isytyp(2,434))
      call getall(Work(1),length,irrepx,434)
      call checksum("List-434:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,435),isytyp(2,435))
      call getall(Work(1),length,irrepx,435)
      call checksum("List-435:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,436),isytyp(2,436))
      call getall(Work(1),length,irrepx,436)
      call checksum("List-436:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,437),isytyp(2,437))
      call getall(Work(1),length,irrepx,437)
      call checksum("List-437:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,438),isytyp(2,438))
      call getall(Work(1),length,irrepx,438)
      call checksum("List-438:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,439),isytyp(2,439))
      call getall(Work(1),length,irrepx,439)
      call checksum("List-439:",Work(1),length)
      Else
      Length=idsymsz(Irrepx,isytyp(1,437),isytyp(2,437))
      call getall(Work(1),length,irrepx,437)
      call checksum("List-437:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,439),isytyp(2,439))
      call getall(Work(1),length,irrepx,439)
      call checksum("List-439:",Work(1),length)
      Endif
#endif 
C
#ifdef _NOSKIP
      If (iuhf .ne.0) then
      Length=idsymsz(Irrepx,isytyp(1,440),isytyp(2,440))
      call getall(Work(1),length,irrepx,440)
      call checksum("List-440:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,441),isytyp(2,441))
      call getall(Work(1),length,irrepx,441)
      call checksum("List-441:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,442),isytyp(2,442))
      call getall(Work(1),length,irrepx,442)
      call checksum("List-442:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,443),isytyp(2,443))
      call getall(Work(1),length,irrepx,443)
      call checksum("List-443:",Work(1),length)
      Else
      Length=idsymsz(Irrepx,isytyp(1,441),isytyp(2,441))
      call getall(Work(1),length,irrepx,441)
      call checksum("List-441:",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,443),isytyp(2,443))
      call getall(Work(1),length,irrepx,443)
      call checksum("List-443:",Work(1),length)
      Endif 
#endif 
      elseif (iside .eq. 2) then

      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,461),isytyp(2,461))
      call getall(Work(1),length,irrepx,461)
      call checksum("List-461;Z(AB,IJ):",Work(1),length)
      Length=idsymsz(Irrepx,isytyp(1,462),isytyp(2,462))
      call getall(Work(1),length,irrepx,462462462462)
      call checksum("List-462;Z(ab,ij):",Work(1),length)
      Endif
      Length=idsymsz(Irrepx,isytyp(1,463),isytyp(2,463))
      call getall(Work(1),length,irrepx,463)
      call checksum("List-463;Z(Ai,Bj):",Work(1),length)
      Endif 
      Return 
      End
