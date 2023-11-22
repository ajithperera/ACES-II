      Subroutine Check_t(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"

      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)

C First T1(A,I)

      IRREPX =1 
      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)
      CALL Getlst(Work,1,1,1,1,90)
      Call Checksum("T1AI", Work, HPA_LENGTH,S)
      IF (UHF) Call Getlst(Work,1,1,1,2,90)
      IF (UHF) Call Checksum("T1ai", Work, HPB_LENGTH,S)

C T2(IJ,AB)
  
      Write(*,*) 
      IF (UHF) Then
         AAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
         BBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         Call Getall(Work, AAAA_LENGTH_IJAB, IRREPX, 44)
         Call Checksum("T2AAAA_ABIJ", Work, AAAA_LENGTH_IJAB,S)
         Call Getall(Work, BBBB_LENGTH_IJAB, IRREPX, 45)
         Call Checksum("T2BBBB_abij", Work, BBBB_LENGTH_IJAB,S)
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 46)
         Call Checksum("T2ABAB_AbIj", Work, ABAB_LENGTH_IJAB,S)
      Else
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 46)
         Call Checksum("T2ABAB_AbIj", Work, ABAB_LENGTH_IJAB,S)
      Endif

      Write(6,*)
      Return
      End


