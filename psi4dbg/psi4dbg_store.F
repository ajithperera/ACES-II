      Subroutine Psi4dbg_store(Work,Length,Iuhf)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)

      Logical null 
      Data onem,One /-1.0,1.0/
#include "sympop.com"
#include "syminf.com"
#include "parcc.com"

C T2-ABIJ
      Irrepx = 1
      Length_44=IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
      If (Iuhf .NE. 0) Length_45=IDSYMSZ(IRREPX,ISYTYP(1,45),
     +                                   ISYTYP(2,45))
      Length_46=IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))

      Write(6,*) 
      Call Getall(Work, Length_44, Irrepx, 44)
      Call Putall(Work, Length_44, Irrepx, 61)
      Call checksum("T2-AAAA:",Work,Length_44)

      If (Iuhf .ne. 0) Then
      Call Getall(Work, Length_45, Irrepx, 45)
      Call Putall(Work, Length_45, Irrepx, 62)
      Call checksum("T2-BBBB:",Work,Length_45)
      Endif 
      Call Getall(Work, Length_46, Irrepx, 46)
      Call Putall(Work, Length_46, Irrepx, 63)
      Call checksum("T2-ABAB:",Work,Length_46)

      Write(6,*) 

      Return
      End

