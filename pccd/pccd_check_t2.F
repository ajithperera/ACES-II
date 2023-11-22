      Subroutine Pccd_check_T2(Work,Length,Iuhf)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)

      Logical null 
      Data onem,One /-1.0,1.0/
#include "sympop.com"
#include "syminf.com"
#include "parcc.com"

C T2-ABIJ
      Irrepx = 1
      Length_63=IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))

      Call Getall(Work, Length_63, Irrepx, 46)
      Call checksum("T2-ABAB :",Work,Length_63)

C T1-AI
      Length_aa =Irpdpd(Irrepx,9)
      call getlst(Work,1,1,1,1,90)
      Call checksum("T1-AA   :",Work,Length_aa)
      If (Iuhf .Ne. 0) Then
      Length_bb =Irpdpd(Irrepx,10)
      call getlst(Work,1,1,1,2,90)
      Call checksum("T1-BB   :",Work,Length_bb)
      Endif
      Write(6,*)

      Return
      End

