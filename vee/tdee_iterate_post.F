      Subroutine tdee_iterate_post(work,Memleft,Iuhf,Irrepx,Iside)

      Implicit Integer(A-Z)

      Dimension Work(Memleft)
      Logical Init 

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"
   
C Source and target singels

      pha_length = Irpdpd(Irrepx,9)
      phb_length = Irpdpd(Irrepx,10)
      hha_length = Irpdpd(Irrepx,21)
      hhb_length = Irpdpd(Irrepx,22)
      ppa_length = Irpdpd(Irrepx,19)
      ppb_length = Irpdpd(Irrepx,20)
    
      Write(6,*) 
      Write(6,"(a,l)") " Entered tdee_iterate-post:"
      Write(6,*) 

      If (Iside .EQ.1) Then
         Call Getlst(Work,1,1,1,1,394)
         Call checksum("list-394-AA      :",Work,Pha_length,s)
         Call Putlst(Work,1,1,1,3,490)
         If (Iuhf .EQ. 1) Then
            Call Getlst(Work,1,1,1,2,394)
            Call checksum("list-394-BB      :",Work,Phb_length,s)
            Call Putlst(Work,1,1,1,4,490)
         Endif 
      Else if (Iside .EQ. 2) Then
         Call Getlst(Work,1,1,1,1,396)
         Call checksum("list-396-AA      :",Work,Pha_length,s)
         Call Putlst(Work,1,1,1,3,490)
         If (Iuhf .EQ. 1) Then
            Call Getlst(Work,1,1,1,2,396)
            Call checksum("list-396-BB      :",Work,Phb_length,s)
            Call Putlst(Work,1,1,1,4,490)
         Endif
      Endif 

      If (Iside .EQ. 1) Then
          If (Iuhf .EQ. 0) Then
              Length = idsymsz(Irrepx,isytyp(1,336),isytyp(2,336))
              Call Getall(Work,Length,Irrepx,336)
              Call checksum("List-336;Z(Ai,Bj):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,463)
           Else 
              Length = idsymsz(Irrepx,isytyp(1,334),isytyp(2,334))
              Call Getall(Work,Length,Irrepx,334)
              Call checksum("List-334;Z(AB,IJ):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,461)
              Length = idsymsz(Irrepx,isytyp(1,335),isytyp(2,335))
              Call Getall(Work,Length,Irrepx,335)
              Call checksum("List-335;Z(ab,ij):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,462)
              Length = idsymsz(Irrepx,isytyp(1,336),isytyp(2,336))
              Call Getall(Work,Length,Irrepx,336)
              Call checksum("List-336;Z(Ai,Bj):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,463)
           Endif
      Else if (Iside .EQ.2) Then
           If (Iuhf .EQ. 0) Then
              Length = idsymsz(Irrepx,isytyp(1,346),isytyp(2,346))
              Call Getall(Work,Length,Irrepx,346)
              Call checksum("List-346;Z(Ai,Bj):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,463)
           Else
              Length = idsymsz(Irrepx,isytyp(1,344),isytyp(2,344))
              Call Getall(Work,Length,Irrepx,344)
              Call checksum("List-344;Z(AB,IJ):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,461)
              Length = idsymsz(Irrepx,isytyp(1,345),isytyp(2,345))
              Call Getall(Work,Length,Irrepx,345)
              Call checksum("List-345;Z(ab,ij):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,462)
              Length = idsymsz(Irrepx,isytyp(1,346),isytyp(2,346))
              Call Getall(Work,Length,Irrepx,346)
              Call checksum("List-346;Z(Ai,Bj):",Work,Length,s)
              Call Putall(Work,Length,Irrepx,463)
           Endif

      Endif 
 
      Return
      End
 
      
