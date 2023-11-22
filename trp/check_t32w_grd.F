      Subroutine Check_t32w_grad(Work,Length,Iuhf)
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
      Length_7=IDSYMSZ(IRREPX,ISYTYP(1,162),ISYTYP(2,162))
      Length_8=IDSYMSZ(IRREPX,ISYTYP(1,163),ISYTYP(2,163))
      Length_9=IDSYMSZ(IRREPX,ISYTYP(1,164),ISYTYP(2,164))
      Length_10=IDSYMSZ(IRREPX,ISYTYP(1,165),ISYTYP(2,165))

      Call Getall(Work, Length_7, Irrepx, 162)
      Call checksum("AAAA_IJKA",Work,Length_7)
      Call Getall(Work, Length_8, Irrepx, 163)
      Call checksum("BBBB_IJKA",Work,Length_8)
      Call Getall(Work, Length_9, Irrepx, 164)
      Call checksum("ABAB-IJKA",Work,Length_9)
      Call Getall(Work, Length_10, Irrepx, 165)
      Call checksum("ABAB-IJAK",Work,Length_10)
      Else
CSSS      Length_7=IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
      Length_10=IDSYMSZ(IRREPX,ISYTYP(1,165),ISYTYP(2,165))
CSSS      Call Getall(Work, Length_7, Irrepx, 107)
CSSS      Call checksum("AAAA_IJKA",Work,Length_7)
      Call Getall(Work, Length_10, Irrepx, 165)
      Call checksum("ABAB-IJAK",Work,Length_10)
      Endif
C ABCI
      If (Uhf) then
      Length_127=IDSYMSZ(IRREPX,ISYTYP(1,166),ISYTYP(2,166))
      Length_128=IDSYMSZ(IRREPX,ISYTYP(1,167),ISYTYP(2,167))
      Length_129=IDSYMSZ(IRREPX,ISYTYP(1,168),ISYTYP(2,168))
      Length_130=IDSYMSZ(IRREPX,ISYTYP(1,169),ISYTYP(2,169))

      Call Getall(Work, Length_127, Irrepx, 166)
      Call checksum("AAAA_ABCI",Work,Length_127)
      Call Getall(Work, Length_128, Irrepx, 167)
      Call checksum("BBBB_ABCI",Work,Length_128)
      Call Getall(Work, Length_129, Irrepx, 168)
      Call checksum("ABAB-ABCI",Work,Length_129)
CSSS      Write(6,"(5(1x,F15.9))") (Work(i),i=1,Length_129)
      Call Getall(Work, Length_130, Irrepx, 169)
      Call checksum("ABAB-ABIC",Work,Length_130)
CSSS      Write(6,"(5(1x,F15.9))") (Work(i),i=1,Length_130)

      Else
CSSS      Length_27=IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,168))
      Length_130=IDSYMSZ(IRREPX,ISYTYP(1,169),ISYTYP(2,169))
CSSS      Call Getall(Work, Length_27, Irrepx, 168)
CSSS      Call checksum("AAAA_ABCI",Work,Length_127)
      Call Getall(Work, Length_130, Irrepx, 169)
      Call checksum("ABAB-ABIC",Work,Length_130)
CSSS      Write(6,"(5(1x,F15.9))") (Work(i),i=1,Length_130)
      Endif 

      Return 
      end 
       

