      Subroutine Pccd_Checkintms(Work,Maxcor,Iuhf,flag)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"
#include "flags.com"

      UHF    = (IUhf .EQ. 1)
      Ao_lad = (Iflags(93) .Gt. 0)

C First Hbar(A,E), Hbar(M,I) and Hbar(M,E) 

      IRREPX =1 
      HHA_LENGTH = IRPDPD(IRREPX,21)
      HHB_LENGTH = IRPDPD(IRREPX,22)
      CALL Getlst(Work,1,1,1,1,91)
      Call Checksum("FMI      ", Work, HHA_LENGTH)
      IF (UHF) Call Getlst(Work,1,1,1,2,91)
      IF (UHF) Call Checksum("Fmi      ", Work, HHB_LENGTH)

      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)
      CALL Getlst(Work,1,1,1,1,93)
      Call Checksum("FME      ", Work, HPA_LENGTH)
      IF (UHF) Call Getlst(Work,1,1,1,2,93)
      IF (UHF) Call Checksum("Fme      ", Work, HPB_LENGTH)

      PPA_LENGTH = IRPDPD(IRREPX,19)
      PPB_LENGTH = IRPDPD(IRREPX,20)
      CALL Getlst(Work,1,1,1,1,92)
      Call Checksum("FEA      ", Work, PPA_LENGTH)
      IF (UHF) Call Getlst(Work,1,1,1,2,92)
      IF (UHF) Call Checksum("Fea      ", Work, PPB_LENGTH)

      If (flag .eq.2) then
         CALL Getlst(Work,1,1,1,9,91)
         Call Checksum("1/4FMI      ", Work, HHA_LENGTH)
         IF (UHF) Call Getlst(Work,1,1,1,10,91)
         IF (UHF) Call Checksum("Fmi      ", Work, HHB_LENGTH)

         CALL Getlst(Work,1,1,1,9,92)
         Call Checksum("1/4FEA      ", Work, PPA_LENGTH)
         IF (UHF) Call Getlst(Work,1,1,1,10,92)
         IF (UHF) Call Checksum("1/4Fea      ", Work, PPB_LENGTH)
      Endif 

C Hbar(MN,IJ)
  
      Write(*,*) 
      IF (UHF) Then
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Checksum("ABAB_MNIJ", Work, ABAB_LENGTH_MNIJ)
      Else
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Checksum("ABAB_MNIJ", Work, ABAB_LENGTH_MNIJ)
      Endif

C Hbar(ABCD)
 
      Write(*,*)
      IF (UHF) Then
         If (.Not. Ao_lad) Then
         ABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))
         Call Getall(Work, ABAB_LENGTH_ABCD, IRREPX, 233)
         Call Checksum("ABAB_ABCD", Work, ABAB_LENGTH_ABCD)
         Endif 
      Else
         If (.Not. Ao_lad) Then
         ABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))
         Call Getall(Work, ABAB_LENGTH_ABCD, IRREPX, 233)
         Call Checksum("ABAB_ABCD", Work, ABAB_LENGTH_ABCD)
         Endif 
      Endif

C Hbar(AI,BC)
      Write(*,*) 
      IF (UHF) Then
         ABAB_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,29),ISYTYP(2,29))
         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))
         Call Getall(Work, ABAB_LENGTH_AIBC, IRREPX, 29)
         Call Checksum("ABAB_AIBC", Work, ABAB_LENGTH_AIBC)
         Call Getall(Work, ABBA_LENGTH_AIBC, IRREPX, 30) 
         Call Checksum("ABBA_AIBC", Work, ABBA_LENGTH_AIBC)
      Else
CSSS         AAAA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
CSSS         Call Getall(Work, AAAA_LENGTH_AIBC, IRREPX, 27)
CSSS         Call Checksum("AAAA_AIBC", Work, AAAA_LENGTH_AIBC)
         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))
         Call Getall(Work, ABBA_LENGTH_AIBC, IRREPX, 30) 
         Call Checksum("ABBA_AIBC", Work, ABBA_LENGTH_AIBC)
      Endif 

C Hbar(IJ,KA)
      Write(*,*)
      IF (UHF) Then
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,9),ISYTYP(2,9))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))
         Call Getall(Work, ABAB_LENGTH_AIBC, IRREPX, 9)
         Call Checksum("ABAB_IJKA", Work, ABAB_LENGTH_IJKA)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 10)
         Call Checksum("ABBA_IJKA", Work, ABBA_LENGTH_IJKA)
      Else
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 10)
         Call Checksum("ABBA_IJKA", Work, ABBA_LENGTH_IJKA)
      Endif
C(MB,EJ)
      Write(6,*)
      IF (UHF) Then
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, ABAB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
      Else
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, ABAB_LENGTH_MBEJ)
      Endif
C (IJ,KL)
         Write(6,*)
         IF (UHF) THEN
         AAAA_LENTH_IJKL =   IDSYMSZ(IRREPX,ISYTYP(1,11),ISYTYP(2,11))
         Call Getall(Work, ABAB_LENGTH_IJKL, IRREPX, 11)
         Call Checksum("ABAB_IJKL", Work, ABAB_LENGTH_IJKL)
         BBBB_LENTH_IJKL =   IDSYMSZ(IRREPX,ISYTYP(1,12),ISYTYP(2,12))
         Call Getall(Work, BBBB_LENGTH_IJKL, IRREPX, 12)
         ELSE
         ABAB_LENGTH_IJKL = IDSYMSZ(IRREPX,ISYTYP(1,13),ISYTYP(2,13))
         Call Getall(Work, ABAB_LENGTH_IJKL, IRREPX, 13)
         Call Checksum("ABAB_IJKL", Work, ABAB_LENGTH_IJKL)
         ENDIF 

      Return
      End


