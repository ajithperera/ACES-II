      Subroutine Check_ints(Work,Length,Iuhf,null)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)

      Logical null 
      Data onem,One /-1.0,1.0/
#include "sympop.com"
#include "syminf.com"
#include "parcc.com"

#ifdef _DEBUG_LVL0
      Irrepx = 1
      If (Ispar .and. Coulomb) Then
         nrow=Isytyp(1,119)
         ncol=Isytyp(2,119)
         Nsize = Idsymsz(Irrepx,Isytyp(1,119),Isytyp(2,119))
         I000 = 1
         I010 = I000 + Nsize
         call getall(Work(I000),Nsize,irrepx,119)
         Call Checksum("119 :",Work(I000),Nsize)
CSS         call output(Work(I000),1,nrow,1,ncol,nrow,ncol,1)
         Nsize = Idsymsz(Irrepx,Isytyp(1,18),Isytyp(2,18))
         I020 = I010 + Nsize
          call getall(Work(i010),Nsize,irrepx,18)
CSS         call output(Work(I010),1,nrow,1,ncol,nrow,ncol,1)
          Call Checksum("21 :",Work(I010),Nsize)
          Call Daxpy(Nsize,Onem,Work(I010),1,Work(I000),1)
          Call Checksum("21 vs 119 :",Work(I000),Nsize)
      Endif 

      Length_23=IDSYMSZ(IRREPX,ISYTYP(1,23),ISYTYP(2,23))
      Length_24=IDSYMSZ(IRREPX,ISYTYP(1,24),ISYTYP(2,24))
      Length_25=IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
      Length_26=IDSYMSZ(IRREPX,ISYTYP(1,26),ISYTYP(2,26))
      Length_17=IDSYMSZ(IRREPX,ISYTYP(1,17),ISYTYP(2,17))
      Length_18=IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
      Call Getall(Work, Length_23, Irrepx, 23)
      Call checksum("List-23",Work,Length_23,s)
      Call Getall(Work, Length_24, Irrepx, 24)
      Call checksum("List-24",Work,Length_24,s)
      Call Getall(Work, Length_25, Irrepx, 25)
      Call checksum("List-25",Work,Length_25,s)
      Call Getall(Work, Length_26, Irrepx, 26)
      Call checksum("List-26",Work,Length_26,s)

      Length_7=IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
      Length_8=IDSYMSZ(IRREPX,ISYTYP(1,8),ISYTYP(2,8))
      Length_9=IDSYMSZ(IRREPX,ISYTYP(1,9),ISYTYP(2,9))
      Length_10=IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

      Call Getall(Work, Length_7, Irrepx, 7)
      Call checksum("List-7 ",Work,Length_7,s)
      Call Getall(Work, Length_8, Irrepx, 8)
      Call checksum("List-8 ",Work,Length_8,s)
      Call Getall(Work, Length_9, Irrepx, 9)
      Call checksum("List-9 ",Work,Length_9,s)
      Call Getall(Work, Length_10, Irrepx, 10)
      Call checksum("List-10",Work,Length_10,s)

C      write(6,*)
C      do irrepr = 1, Nirrep
C         irrepl = dirprd(irrepr,irrepx)
C         nrow = irpdpd(irrepl,isytyp(1,10))
C         ncol = irpdpd(irrepr,isytyp(2,10))
C         write(6,"(a,2i4)") "nrow,ncol:",nrow,ncol
C         call getlst(work,1,ncol,1,irrepr,10)
C         call checksum("list-10",work,nrow*ncol,s)
C      enddo
C      write(6,*)

      Length_27=IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
      Length_28=IDSYMSZ(IRREPX,ISYTYP(1,28),ISYTYP(2,28))
      Length_29=IDSYMSZ(IRREPX,ISYTYP(1,29),ISYTYP(2,29))
      Length_30=IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))

      Call Getall(Work, Length_27, Irrepx, 27)
      Call checksum("List-27",Work,Length_27,s)
      Call Getall(Work, Length_28, Irrepx, 28)
      Call checksum("List-28",Work,Length_28,s)
      Call Getall(Work, Length_29, Irrepx, 29)
      Call checksum("List-29",Work,Length_29,s)
      Call Getall(Work, Length_30, Irrepx, 30)
      Call checksum("List-30",Work,Length_30,s)
#ifdef _NOSKIP
      if (null) then 
         Call getall(Work, Length_17, Irrepx, 17)
         call checksum("17",work,Length_17)
         Call Putall(Work, Length_17, Irrepx, 147)
         Call getall(Work, Length_18, Irrepx, 18)
         call checksum("18",work,Length_18)
         Call Putall(Work, Length_18, Irrepx, 148)

         Write(6,*) "Zeroing 17 and 18"
         call dzero(Work,length_17)
         Call Putall(Work, Length_17, Irrepx, 121)
         call dzero(Work,length_18)
         Call Putall(Work, Length_18, Irrepx, 122)
         call checksum("17",work,Length_17)
         call checksum("18",work,Length_18)
      else 
         Write(6,*) "Restoring 17 and 18"
         Call getall(Work, Length_17, Irrepx, 147)
         call checksum("17",work,Length_17)
         Call Putall(Work, Length_17, Irrepx, 121)
         Call Getall(Work, Length_18, Irrepx, 148)
         call checksum("18",work,Length_18)
         Call Putall(Work, Length_18, Irrepx, 122)
      endif 
#endif
       
#endif
      Return 
      end 
       

