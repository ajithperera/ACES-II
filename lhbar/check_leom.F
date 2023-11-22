      Subroutine Check_leom(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"

      UHF = (IUhf .EQ. 1)

C First L(ia)) 

      IRREPX =1 
      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)
      CALL Getlst(Work,1,1,1,3,90)
      Call Checksum("AA_L1    :", Work, HPA_LENGTH)
      IF (UHF) Call Getlst(Work,1,1,1,4,90)
      IF (UHF) Call Checksum("BB-L1    :", Work, HPB_LENGTH)

C L(ij,ab) increment 

      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,61),ISYTYP(2,61))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,62),ISYTYP(2,62))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))
         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 61)
         Call Checksum("AAAA_L2  :", Work, AAAA_LENGTH_IJKA)
         Call Getall(Work, BBBB_LENGTH_IJKA, IRREPX, 62)
         Call Checksum("BBBB_L2  :", Work, BBBB_LENGTH_IJKA)
         Call Getall(Work, ABAB_LENGTH_IJKA, IRREPX, 63)
         Call Checksum("ABAB_L2  :", Work, ABAB_LENGTH_IJKA)
      Else
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX,63) 
         Call Checksum("ABBA_L2  :", Work, ABBA_LENGTH_IJKA)
      Endif

C L(ij,ab) 

      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,144),ISYTYP(2,144))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,145),ISYTYP(2,145))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,146),ISYTYP(2,146))
         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 144)
         Call Checksum("AAAA_L2  :", Work, AAAA_LENGTH_IJKA)
         Call Getall(Work, BBBB_LENGTH_IJKA, IRREPX, 145)
         Call Checksum("BBBB_L2  :", Work, BBBB_LENGTH_IJKA)
         Call Getall(Work, BBBB_LENGTH_IJKA, IRREPX, 145)
         Call Getall(Work, ABAB_LENGTH_IJKA, IRREPX, 146)
         Call Checksum("ABAB_L2  :", Work, ABAB_LENGTH_IJKA)
      Else
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,146),ISYTYP(2,146))
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX,146)
         Call Checksum("ABBA_L2  :", Work, ABBA_LENGTH_IJKA)
      Endif


      Return
      End
