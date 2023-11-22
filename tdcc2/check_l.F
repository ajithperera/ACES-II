      Subroutine Check_l(Work,Maxcor,Iuhf)


      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"

      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)

C First L1(A,I)

      IRREPX =1 
      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)
      CALL Getlst(Work,1,1,1,1,190)
      Call Checksum("L1IA", Work, HPA_LENGTH,S)
      IF (UHF) Call Getlst(Work,1,1,1,2,190)
      IF (UHF) Call Checksum("L1ia", Work, HPB_LENGTH,S)

C L2(IJ,AB)
  
      Write(*,*) 
      IF (UHF) Then
         AAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,144),ISYTYP(2,144))
         BBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,145),ISYTYP(2,145))
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,146),ISYTYP(2,146))
         Call Getall(Work, AAAA_LENGTH_IJAB, IRREPX, 144)
         Call Checksum("L2AAAA_IJAB", Work, AAAA_LENGTH_IJAB,S)
         Call Getall(Work, BBBB_LENGTH_IJAB, IRREPX, 145)
         Call Checksum("L2BBBB_ijab", Work, BBBB_LENGTH_IJAB,S)
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 146)
         Call Checksum("L2ABAB_IjAb", Work, ABAB_LENGTH_IJAB,S)
      Else
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,146),ISYTYP(2,146))
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 146)
         Call Checksum("L2ABAB_IjAb", Work, ABAB_LENGTH_IJAB,S)
      Endif

      Write(6,*)

      Return
      End


