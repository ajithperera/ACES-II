      Subroutine Check_vee_lists_init(Work,Maxcor,Irrepx,Iuhf,Iside)
      Implicit Integer(A-Z)

      Double Precision Work(Maxcor)
      Logical UHF
      Logical David,Block_david,Lanczos

#include "machsp.com"
#include "sympop.com"
      Common /Rootsearch/David,Block_david,Lanczos

      pha_length = Irpdpd(Irrepx,9)
      phb_length = Irpdpd(Irrepx,10)
      hha_length = Irpdpd(Irrepx,21)
      hhb_length = Irpdpd(Irrepx,22)
      ppa_length = Irpdpd(Irrepx,19)
      ppb_length = Irpdpd(Irrepx,20)

      Write(6,*) 
      If (iside .Eq. 1) then 

      Call Getlst(Work,1,1,1,1,490)
      Call checksum("list-490-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,2,490)
         Call checksum("list-490-BB      :",Work,Phb_length,s)
      Endif

      Elseif (Iside .Eq. 2) Then

      If (Lanczos) then
      Call Getlst(Work,1,1,1,1,493)
      Call checksum("list-493-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,2,493)
         Call checksum("list-493-BB      :",Work,Phb_length,s)
      Endif
      Call Getlst(Work,1,1,1,1,490)
      Call checksum("list-490-AA      :",Work,Pha_length,s)
      If (Iuhf .Ne. 0) Then
         Call Getlst(Work,1,1,1,2,490)
         Call checksum("list-490-BB      :",Work,Phb_length,s)
      Endif 
      Endif 
      Endif 

      If (iside .eq. 1) Then 
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,444),isytyp(2,444))
      call getall(Work(1),length,irrepx,444)
      call checksum("List-444;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,445),isytyp(2,445))
      call getall(Work(1),length,irrepx,445)
      call checksum("List-445;Z(ab,ij):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
      call getall(Work(1),length,irrepx,446)
      call checksum("List-446;Z(ab,ij):",Work(1),length,s)
      Else
      Length=idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
      call getall(Work(1),length,irrepx,446)
      call checksum("List-446;Z(Ai,Bj):",Work(1),length,s)
      Endif 

      elseif (iside .eq. 2) then

      If(Lanczos) then
      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,448),isytyp(2,448))
      call getall(Work(1),length,irrepx,448)
      call checksum("List-448;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,449),isytyp(2,449))
      call getall(Work(1),length,irrepx,449)
      call checksum("List-449;Z(ab,ij):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,450),isytyp(2,450))
      call getall(Work(1),length,irrepx,450)
      call checksum("List-450;Z(ab,ij):",Work(1),length,s)
      Else 
      Length=idsymsz(Irrepx,isytyp(1,450),isytyp(2,450))
      call getall(Work(1),length,irrepx,450)
      call checksum("List-450;Z(Ai,Bj):",Work(1),length,s)
      Endif
      Else

      If (iuhf .ne. 0) Then
      Length=idsymsz(Irrepx,isytyp(1,444),isytyp(2,444))
      call getall(Work(1),length,irrepx,444)
      call checksum("List-444;Z(AB,IJ):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,445),isytyp(2,445))
      call getall(Work(1),length,irrepx,445)
      call checksum("List-445;Z(ab,ij):",Work(1),length,s)
      Length=idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
      call getall(Work(1),length,irrepx,446)
      call checksum("List-446;Z(ab,ij):",Work(1),length,s)
      Else
      Length=idsymsz(Irrepx,isytyp(1,446),isytyp(2,446))
      call getall(Work(1),length,irrepx,446)
      call checksum("List-446;Z(Ai,Bj):",Work(1),length,s)
      Endif 
      Endif 
      Endif 

      Return 
      End
