      Subroutine Check_fintms(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"

      UHF = (IUhf .EQ. 1)

C First Hbar(A,E), Hbar(M,I) and Hbar(M,E) 

      IRREPX =1 
      HHA_LENGTH = IRPDPD(IRREPX,21)
      HHB_LENGTH = IRPDPD(IRREPX,22)
      CALL Getlst(Work,1,1,1,1,91)
      Call Checksum("FMI", Work, HHA_LENGTH)
      CALL Getlst(Work,1,1,1,9,91)
      Call Checksum("FMI", Work, HHA_LENGTH)
      IF (UHF) Call Getlst(Work,1,1,1,2,91)
      IF (UHF) Call Checksum("Fmi", Work, HHB_LENGTH)
      IF (UHF) Call Getlst(Work,1,1,1,10,91)
      IF (UHF) Call Checksum("Fmi", Work, HHB_LENGTH)

      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)
      CALL Getlst(Work,1,1,1,1,93)
      Call Checksum("FME", Work, HPA_LENGTH)
      IF (UHF) Call Getlst(Work,1,1,1,2,93)
      IF (UHF) Call Checksum("Fme", Work, HPB_LENGTH)

      PPA_LENGTH = IRPDPD(IRREPX,19)
      PPB_LENGTH = IRPDPD(IRREPX,20)
      CALL Getlst(Work,1,1,1,1,92)
      Call Checksum("FEA", Work, PPA_LENGTH)
      CALL Getlst(Work,1,1,1,9,92)
      Call Checksum("FEA", Work, PPA_LENGTH)
      IF (UHF) Call Getlst(Work,1,1,1,2,92)
      IF (UHF) Call Checksum("Fea", Work, PPB_LENGTH)
      IF (UHF) Call Getlst(Work,1,1,1,10,92)
      IF (UHF) Call Checksum("Fea", Work, PPB_LENGTH)

      Return
      End


