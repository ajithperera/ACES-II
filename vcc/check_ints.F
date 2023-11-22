      Subroutine Check_ints(Work,Length,Iuhf,null)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)

      Logical null,uhf, Ao_lad,noabcd
      Data onem,One /-1.0,1.0/
#include "sympop.com"
#include "syminf.com"
#include "parcc.com"
#include "flags.com"

      UHF    = (IUhf .EQ. 1)
      AO_lad = (Iflags(93) .Gt. 0) 
      Noabcd = (Iflags(1) .eq. 1)

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
C IJKL
      If (Uhf) Then
      Length_11=IDSYMSZ(IRREPX,ISYTYP(1,11),ISYTYP(2,11))
      Length_12=IDSYMSZ(IRREPX,ISYTYP(1,12),ISYTYP(2,12))
      Length_13=IDSYMSZ(IRREPX,ISYTYP(1,13),ISYTYP(2,13))

      Call Getall(Work, Length_11, Irrepx, 11)
      Call checksum("IJKL-AAAA",Work,Length_11)
CSSS      Write(6,"(6(1x,F15.10))") (Work(i),i=1,Length_11)
      Call Getall(Work, Length_12, Irrepx, 12)
      Call checksum("IJKL-BBBB",Work,Length_12)
      Call Getall(Work, Length_13, Irrepx, 13)
CSSS      Write(6,"(6(1x,F15.10))") (Work(i),i=1,Length_13)
      Call checksum("IJKL-ABAB",Work,Length_13)
      Else
      Length_13=IDSYMSZ(IRREPX,ISYTYP(1,13),ISYTYP(2,13))
      Call Getall(Work, Length_13, Irrepx, 13)
      Call checksum("IJKL-ABAB",Work,Length_13)
      Endif 
C ABIJ
      If (uhf) Then
      Length_14=IDSYMSZ(IRREPX,ISYTYP(1,14),ISYTYP(2,14))
      Length_15=IDSYMSZ(IRREPX,ISYTYP(1,15),ISYTYP(2,15))
      Length_16=IDSYMSZ(IRREPX,ISYTYP(1,16),ISYTYP(2,16))
      Length_17=IDSYMSZ(IRREPX,ISYTYP(1,17),ISYTYP(2,17))
      Length_18=IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
      Length_19=IDSYMSZ(IRREPX,ISYTYP(1,19),ISYTYP(2,19))
      Length_20=IDSYMSZ(IRREPX,ISYTYP(1,20),ISYTYP(2,20))
      Length_21=IDSYMSZ(IRREPX,ISYTYP(1,21),ISYTYP(2,21))
      Length_22=IDSYMSZ(IRREPX,ISYTYP(1,22),ISYTYP(2,22))

      Call Getall(Work, Length_14, Irrepx, 14)
      Call checksum("ABIJ-AAAA",Work,Length_14)
      Call Getall(Work, Length_15, Irrepx, 15)
      Call checksum("ABIJ-BBBB",Work,Length_15)
      Call Getall(Work, Length_16, Irrepx, 16)
      Call checksum("ABIJ-ABAB",Work,Length_16)
      Call Getall(Work, Length_17, Irrepx, 17)
      Call checksum("ABIJ-BBAA",Work,Length_17)
      Call Getall(Work, Length_18, Irrepx, 18)
      Call checksum("ABIJ-AABB",Work,Length_18)
      Call Getall(Work, Length_19, Irrepx, 19)
      Call checksum("ABIJ-AAAA",Work,Length_19)
      Call Getall(Work, Length_20, Irrepx, 20)
      Call checksum("ABIJ-BBBB",Work,Length_20)
      Call Getall(Work, Length_21, Irrepx, 21)
      Call checksum("ABIJ-ABBA",Work,Length_21)
      Call Getall(Work, Length_22, Irrepx, 22)
      Call checksum("ABIJ-BAAB",Work,Length_22)
      Else
      Length_16=IDSYMSZ(IRREPX,ISYTYP(1,16),ISYTYP(2,16))
      Length_18=IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
      Length_19=IDSYMSZ(IRREPX,ISYTYP(1,19),ISYTYP(2,19))
      Length_21=IDSYMSZ(IRREPX,ISYTYP(1,21),ISYTYP(2,21))
      
      Call Getall(Work, Length_16, Irrepx, 16)
      Call checksum("ABIJ-ABAB",Work,Length_16)
      Call Getall(Work, Length_18, Irrepx, 18)
      Call checksum("ABIJ-AABB",Work,Length_18)
      Call Getall(Work, Length_19, Irrepx, 19)
      Call checksum("ABIJ-AAAA",Work,Length_19)
      Call Getall(Work, Length_21, Irrepx, 21)
      Call checksum("ABIJ-ABBA",Work,Length_21)
      Endif 
CMBEJ
      If (uhf) Then
      Length_23=IDSYMSZ(IRREPX,ISYTYP(1,23),ISYTYP(2,23))
      Length_24=IDSYMSZ(IRREPX,ISYTYP(1,24),ISYTYP(2,24))
      Length_25=IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
      Length_26=IDSYMSZ(IRREPX,ISYTYP(1,26),ISYTYP(2,26))

      Call Getall(Work, Length_23, Irrepx, 23)
CSSS      Write(6,"(6(1x,F15.10))") (Work(i),i=1,Length_23)
      Call checksum("MBEJ-AAAA",Work,Length_23)
      Call Getall(Work, Length_24, Irrepx, 24)
      Call checksum("MBEJ-BBBB",Work,Length_24)
      Call Getall(Work, Length_25, Irrepx, 25)
      Call checksum("MBEJ-ABBA",Work,Length_25)
CSSS      Write(6,"(6(1x,F15.10))") (Work(i),i=1,Length_25)
      Call Getall(Work, Length_26, Irrepx, 26)
      Call checksum("MBEJ-BAAB",Work,Length_26)
CSSS      Write(6,"(6(1x,F15.10))") (Work(i),i=1,Length_26)

      If (.Not. AO_lad) Then
      LAAAA_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,231),ISYTYP(2,231))
      LBBBB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,232),ISYTYP(2,232))
      LABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))
        
      Call Getall(Work, LAAAA_LENGTH_ABCD, IRREPX, 231)
      Call Checksum("AAAA_ABCD", Work, LAAAA_LENGTH_ABCD)
      Call Getall(Work, LBBBB_LENGTH_ABCD, IRREPX, 232)
      Call Checksum("BBBB_ABCD", Work, LBBBB_LENGTH_ABCD)
      Call Getall(Work, LABAB_LENGTH_ABCD, IRREPX, 233)
      Call Checksum("ABAB_ABCD", Work, LABAB_LENGTH_ABCD)
      Endif 

      Else

      Length_23=IDSYMSZ(IRREPX,ISYTYP(1,23),ISYTYP(2,23))
      Length_25=IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))

      If (.Not. AO_lad) Then
      LABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))
      Endif 

      Call Getall(Work, Length_23, Irrepx, 23)
      Call checksum("MBEJ-AAAA",Work,Length_23)
      Call Getall(Work, Length_25, Irrepx, 25)
      Call checksum("MBEJ-ABBA",Work,Length_25)
     
      If (.Not. AO_lad) Then
      Call Getall(Work, LABAB_LENGTH_ABCD, IRREPX, 233)
      Call Checksum("ABAB_ABCD", Work, LABAB_LENGTH_ABCD)
      Endif 

      Endif 
C IJKA 
      If (Uhf) then
      Length_7=IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
      Length_8=IDSYMSZ(IRREPX,ISYTYP(1,8),ISYTYP(2,8))
      Length_9=IDSYMSZ(IRREPX,ISYTYP(1,9),ISYTYP(2,9))
      Length_10=IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

      Call Getall(Work, Length_7, Irrepx, 7)
      Call checksum("AAAA_IJKA",Work,Length_7)
      Call Getall(Work, Length_8, Irrepx, 8)
      Call checksum("BBBB_IJKA",Work,Length_8)
      Call Getall(Work, Length_9, Irrepx, 9)
      Call checksum("ABAB-IJKA",Work,Length_9)
      Call Getall(Work, Length_10, Irrepx, 10)
      Call checksum("ABAB-IJAK",Work,Length_10)
      Else 
      Length_7=IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
      Length_10=IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))
      Call Getall(Work, Length_7, Irrepx, 7)
      Call checksum("AAAA_IJKA",Work,Length_7)
      Call Getall(Work, Length_10, Irrepx, 10)
      Call checksum("ABAB-IJAK",Work,Length_10)
      Endif 

C ABCI
      If (Uhf) then
      Length_27=IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
      Length_28=IDSYMSZ(IRREPX,ISYTYP(1,28),ISYTYP(2,28))
      Length_29=IDSYMSZ(IRREPX,ISYTYP(1,29),ISYTYP(2,29))
      Length_30=IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))

      Call Getall(Work, Length_27, Irrepx, 27)
      Call checksum("AAAA_ABCI",Work,Length_27)
CSSS      Write(6,"(6(1x,F15.10))") (Work(i),i=1,Length_27)
      Call Getall(Work, Length_28, Irrepx, 28)
      Call checksum("BBBB_ABCI",Work,Length_28)
      Call Getall(Work, Length_29, Irrepx, 29)
      Call checksum("ABAB-ABCI",Work,Length_29)
      Call Getall(Work, Length_30, Irrepx, 30)
      Call checksum("ABAB-ABIC",Work,Length_30)
CSSS      Write(6,"(6(1x,F15.10))") (Work(i),i=1,Length_30)

      Else
CSSS      Length_27=IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
      Length_30=IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))
CSSS      Call Getall(Work, Length_27, Irrepx, 27)
CSSS      Call checksum("AAAA_ABCI",Work,Length_27)
      Call Getall(Work, Length_30, Irrepx, 30)
      Call checksum("ABAB-ABIC",Work,Length_30)
      Endif 

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
       

