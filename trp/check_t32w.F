      Subroutine Check_t32w(Work,Length,Iuhf)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)

      Logical null,uhf 
      Data onem,One /-1.0,1.0/
#include "sympop.com"
#include "syminf.com"
#include "parcc.com"

      UHF = (IUhf .EQ. 1)
      Irrepx = 1
      Write(6,*)
C IJKA 
      If (Uhf) then
      Length_7=IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
      Length_8=IDSYMSZ(IRREPX,ISYTYP(1,8),ISYTYP(2,8))
      Length_9=IDSYMSZ(IRREPX,ISYTYP(1,9),ISYTYP(2,9))
      Length_10=IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

      Call Getall(Work, Length_7, Irrepx, 107)
      Call checksum("AAAA_IJKA",Work,Length_7)
      Call Getall(Work, Length_8, Irrepx, 108)
      Call checksum("BBBB_IJKA",Work,Length_8)
      Call Getall(Work, Length_9, Irrepx, 109)
      Call checksum("ABAB-IJKA",Work,Length_9)
      Call Getall(Work, Length_10, Irrepx, 110)
      Call checksum("ABAB-IJAK",Work,Length_10)
      Else 
      Length_7=IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
      Length_10=IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))
      Call Getall(Work, Length_7, Irrepx, 107)
      Call checksum("AAAA_IJKA",Work,Length_7)
      Call Getall(Work, Length_10, Irrepx, 110)
      Call checksum("ABAB-IJAK",Work,Length_10)
      Endif 

C ABCI
      If (Uhf) then
      Length_27=IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
      Length_28=IDSYMSZ(IRREPX,ISYTYP(1,28),ISYTYP(2,28))
      Length_29=IDSYMSZ(IRREPX,ISYTYP(1,29),ISYTYP(2,29))
      Length_30=IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))

      Call Getall(Work, Length_27, Irrepx, 127)
      Call checksum("AAAA_ABCI",Work,Length_27)
      Call Getall(Work, Length_28, Irrepx, 128)
      Call checksum("BBBB_ABCI",Work,Length_28)
      Call Getall(Work, Length_29, Irrepx, 129)
C      write(6,"(a)") "List-29-UHF"
C      write(6,"(4(1x,F15.9))") (Work(i),i=1,Length_29)
      Call checksum("ABAB-ABCI",Work,Length_29)
      Call Getall(Work, Length_30, Irrepx, 130)
C      write(6,"(a)") "List-30-UHF"
C      write(6,"(4(1x,F15.9))") (Work(i),i=1,Length_29)
      Call checksum("ABAB-ABIC",Work,Length_30)

      Else
CSSS      Length_27=IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
      Length_30=IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))
CSSS      Call Getall(Work, Length_27, Irrepx, 127)
CSSS      Call checksum("AAAA_ABCI",Work,Length_27)
      Call Getall(Work, Length_30, Irrepx, 130)
C      write(6,"(a)") "List-30-RHF"
C      write(6,"(4(1x,F15.9))") (Work(i),i=1,Length_30)
      Call checksum("ABAB-ABIC",Work,Length_30)
      Endif 

      Return 
      end 
       

