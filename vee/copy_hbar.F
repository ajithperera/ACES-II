      Subroutine Copy_Hbar(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"

      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)
      IOFF= 200
      IRREPX = 1

C First Hbar(A,E), Hbar(M,I) and Hbar(M,E) 

      HHA_LENGTH = IRPDPD(IRREPX,21)
      HHB_LENGTH = IRPDPD(IRREPX,22)

      CALL Getlst(Work,1,1,1,1,91+IOFF)
      Call Putlst(Work,1,1,1,1,91)

      If (UHF) Then
          Call Getlst(Work,1,1,1,2,91+IOFF)
          Call putlst(Work,1,1,1,2,91)
      Endif

      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)

      CALL Getlst(Work,1,1,1,1,93+IOFF)
      CALL Putlst(Work,1,1,1,1,93)
      If (UHF) Then
         Call Getlst(Work,1,1,1,2,93+IOFF)
         Call Putlst(Work,1,1,1,2,93) 
      Endif 

      PPA_LENGTH = IRPDPD(IRREPX,19)
      PPB_LENGTH = IRPDPD(IRREPX,20)

      CALL Getlst(Work,1,1,1,1,92+IOFF)
      CALL Putlst(Work,1,1,1,1,92)
      If (UHF) Then
         Call Getlst(Work,1,1,1,2,92+IOFF)
         Call Putlst(Work,1,1,1,2,92)
      Endif 

C Hbar(MN,IJ)
  
      IF (UHF) Then
         AAAA_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,51),ISYTYP(2,51))
         BBBB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,52),ISYTYP(2,52))
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))

         Call Getall(Work, AAAA_LENGTH_MNIJ, IRREPX, 51+IOFF)
         Call Putall(Work, AAAA_LENGTH_MNIJ, IRREPX, 51)
         Call Getall(Work, BBBB_LENGTH_MNIJ, IRREPX, 52+IOFF)
         Call Putall(Work, BBBB_LENGTH_MNIJ, IRREPX, 52)
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53+IOFF)
         Call Putall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53)
         
      Else
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))

         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53+IOFF)
         Call Putall(Work, AAAA_LENGTH_MNIJ, IRREPX, 53)

      Endif

C Hbar(MB,EJ)

      IF (UHF) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
	 ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54+IOFF)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55+IOFF)
         Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56+IOFF)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57+IOFF)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58+IOFF)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59+IOFF)
         Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)

      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54+IOFF)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56+IOFF)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58+IOFF)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)

      Endif

C Hbar(AB,CI)

      IF (UHF) Then
         AAAA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,127))
         BBBB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,128),ISYTYP(2,128))
         ABAB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,129),ISYTYP(2,129))
         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))

         Call Getall(Work, AAAA_LENGTH_ABCI, IRREPX, 127+IOFF)
         Call Putall(Work, AAAA_LENGTH_ABCI, IRREPX, 127)
         Call Getall(Work, BBBB_LENGTH_ABCI, IRREPX, 128+IOFF)
         Call Putall(Work, BBBB_LENGTH_ABCI, IRREPX, 128)
         Call Getall(Work, ABAB_LENGTH_ABCI, IRREPX, 129+IOFF)
         Call Putall(Work, ABAB_LENGTH_ABCI, IRREPX, 129)
         Call Getall(Work, ABBA_LENGTH_ABCI, IRREPX, 130+IOFF)
         Call Putall(Work, ABBA_LENGTH_ABCI, IRREPX, 130)

      Else
CSSS         AAAA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,127))
         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))
 
CSSS         Call Getall(Work, AAAA_LENGTH_ABCI, IRREPX, 127+IOFF)
CSSS         Call Putall(Work, AAAA_LENGTH_ABCI, IRREPX, 127)
         Call Getall(Work, ABBA_LENGTH_ABCI, IRREPX, 130+IOFF)
         Call Putall(Work, ABBA_LENGTH_ABCI, IRREPX, 130)

      Endif

C Hbar(IA,JK)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,107),ISYTYP(2,107))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,108),ISYTYP(2,108))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,109),ISYTYP(2,109))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))

         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 107+IOFF)
         Call Putall(Work, AAAA_LENGTH_IJKA, IRREPX, 107)
         Call Getall(Work, BBBB_LENGTH_IJKA, IRREPX, 108+IOFF)
         Call Putall(Work, BBBB_LENGTH_IJKA, IRREPX, 108)
         Call Getall(Work, ABAB_LENGTH_IJKA, IRREPX, 109+IOFF)
         Call Putall(Work, ABAB_LENGTH_IJKA, IRREPX, 109)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 110+IOFF)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)

      Else
CSSS         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,107),ISYTYP(2,107))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))

CSSS         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 107+IOFF)
CSSS         Call Putall(Work, AAAA_LENGTH_IJKA, IRREPX, 107)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 110+IOFF)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)
      Endif
#ifdef _DEBUG_LVLM
      Call checkhbar(Work,Maxcor,Iuhf)
#endif 
      Return
      End


