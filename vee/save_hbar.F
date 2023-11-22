      Subroutine Save_Hbar(Work,Maxcor,Iuhf)
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

      Call Updmoi(1,HHA_LENGTH,1,91+IOFF,0,0)
      Call Updmoi(1,HHB_LENGTH,2,91+IOFF,0,0)

      CALL Getlst(Work,1,1,1,1,91)
      Call Putlst(Work,1,1,1,1,91+IOFF)

      If (UHF) Then
         Call Getlst(Work,1,1,1,2,91)
         Call Putlst(Work,1,1,1,2,91+IOFF)
      Endif 

      HPA_LENGTH = IRPDPD(IRREPX,9)
      HPB_LENGTH = IRPDPD(IRREPX,10)

      Call Updmoi(1,HPA_LENGTH,1,93+IOFF,0,0)
      Call Updmoi(1,HPB_LENGTH,2,93+IOFF,0,0)

      CALL Getlst(Work,1,1,1,1,93)
      CALL Putlst(Work,1,1,1,1,93+IOFF)

      IF (UHF) Then
         Call Getlst(Work,1,1,1,2,93)
         Call Putlst(Work,1,1,1,2,93+IOFF)
      Endif 

      PPA_LENGTH = IRPDPD(IRREPX,19)
      PPB_LENGTH = IRPDPD(IRREPX,20)

      Call Updmoi(1,PPA_LENGTH,1,92+IOFF,0,0)
      Call Updmoi(1,PPB_LENGTH,2,92+IOFF,0,0)

      CALL Getlst(Work,1,1,1,1,92)
      CALL Putlst(Work,1,1,1,1,92+IOFF)
      If (UHF) Then
         Call Getlst(Work,1,1,1,2,92)
         Call Putlst(Work,1,1,1,2,92+IOFF)
      Endif 

C Hbar(MN,IJ)
  
      IF (UHF) Then
         AAAA_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,51),ISYTYP(2,51))
         BBBB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,52),ISYTYP(2,52))
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))

         Call Inipck(IRREPX,ISYTYP(1,51),ISYTYP(2,51),51+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,52),ISYTYP(2,52),52+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,53),ISYTYP(2,53),53+IOFF,IMODE,0,1)

         Call Getall(Work, AAAA_LENGTH_MNIJ, IRREPX, 51)
         Call Putall(Work, AAAA_LENGTH_MNIJ, IRREPX, 51+IOFF)
         Call Getall(Work, BBBB_LENGTH_MNIJ, IRREPX, 52)
         Call Putall(Work, BBBB_LENGTH_MNIJ, IRREPX, 52+IOFF)
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Putall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53+IOFF)
         
      Else
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))

         Call Inipck(IRREPX,ISYTYP(1,53),ISYTYP(2,53),53+IOFF,IMODE,0,1)

         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Putall(Work, AAAA_LENGTH_MNIJ, IRREPX, 53+IOFF)

      Endif

C Hbar(AI,BC)

      IF (UHF) Then
         AAAA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
         BBBB_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,28),ISYTYP(2,28))
         ABAB_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,29),ISYTYP(2,29))
         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))

         Call Inipck(IRREPX,ISYTYP(1,27),ISYTYP(2,27),27+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,28),ISYTYP(2,28),28+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,29),ISYTYP(2,29),29+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,30),ISYTYP(2,30),30+IOFF,IMODE,0,1)

         Call Getall(Work, AAAA_LENGTH_AIBC, IRREPX, 27)
         Call Putall(Work, AAAA_LENGTH_AIBC, IRREPX, 27+IOFF)
         Call Getall(Work, BBBB_LENGTH_AIBC, IRREPX, 28)
         Call Putall(Work, BBBB_LENGTH_AIBC, IRREPX, 28+IOFF)
         Call Getall(Work, ABAB_LENGTH_AIBC, IRREPX, 29)
         Call Putall(Work, ABAB_LENGTH_AIBC, IRREPX, 29+IOFF)
         Call Getall(Work, ABBA_LENGTH_AIBC, IRREPX, 30) 
         Call Putall(Work, ABBA_LENGTH_AIBC, IRREPX, 30+IOFF)

      Else

         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))

         Call Inipck(IRREPX,ISYTYP(1,30),ISYTYP(2,30),27+IOFF,IMODE,0,1)

         Call Getall(Work, ABBA_LENGTH_AIBC, IRREPX, 30) 
         Call Putall(Work, ABBA_LENGTH_AIBC, IRREPX, 30+IFF) 

      Endif 

C Hbar(IJ,KA)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,8),ISYTYP(2,8))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,9),ISYTYP(2,9))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

         Call Inipck(IRREPX,ISYTYP(1,7),ISYTYP(2,7),7+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,8),ISYTYP(2,8),8+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,9),ISYTYP(2,9),9+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,10),ISYTYP(2,10),10+IOFF,IMODE,0,1)

         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 7)
         Call Putall(Work, AAAA_LENGTH_IJKA, IRREPX, 7+IOFF)
         Call Getall(Work, BBBB_LENGTH_IJKA, IRREPX, 8)
         Call Putall(Work, BBBB_LENGTH_IJKA, IRREPX, 8+IOFF)
         Call Getall(Work, ABAB_LENGTH_IJKA, IRREPX, 9)
         Call Putall(Work, ABAB_LENGTH_IJKA, IRREPX, 9+IOFF)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 10)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 10+IOFF)

      Else
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

         Call Inipck(IRREPX,ISYTYP(1,10),ISYTYP(2,10),10+IOFF,IMODE,0,1)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 10)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 10+IOFF)
      Endif

C Hbar(MB,EJ)

      IF (UHF) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
	 ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))

         Call Inipck(IRREPX,ISYTYP(1,54),ISYTYP(2,54),54+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,55),ISYTYP(2,55),55+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,56),ISYTYP(2,56),56+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,57),ISYTYP(2,57),57+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,58),ISYTYP(2,58),58+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,59),ISYTYP(2,59),59+IOFF,IMODE,0,1)

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54+IOFF)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55+IOFF)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56+IOFF)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57+IOFF)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58+IOFF)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
         Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 59+IOFF)


      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

         Call Inipck(IRREPX,ISYTYP(1,54),ISYTYP(2,54),54+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,56),ISYTYP(2,56),56+IOFF,IMODE,0,1)
         Call Inipck(IRREPX,ISYTYP(1,58),ISYTYP(2,58),58+IOFF,IMODE,0,1)

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54+IOFF)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56+IOFF)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58+IOFF)

      Endif

C Hbar(AB,CI)

      IF (UHF) Then
         AAAA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,127))
         BBBB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,128),ISYTYP(2,128))
         ABAB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,129),ISYTYP(2,129))
         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))

         Call Inipck(IRREPX,ISYTYP(1,127),ISYTYP(2,127),127+IOFF,IMODE,
     +               0,1)
         Call Inipck(IRREPX,ISYTYP(1,128),ISYTYP(2,128),128+IOFF,IMODE,
     +               0,1)
         Call Inipck(IRREPX,ISYTYP(1,129),ISYTYP(2,129),129+IOFF,IMODE,
     +               0,1)
         Call Inipck(IRREPX,ISYTYP(1,130),ISYTYP(2,130),130+IOFF,IMODE,
     +               0,1)

         Call Getall(Work, AAAA_LENGTH_ABCI, IRREPX, 127)
         Call Putall(Work, AAAA_LENGTH_ABCI, IRREPX, 127+IOFF)
         Call Getall(Work, BBBB_LENGTH_ABCI, IRREPX, 128)
         Call Putall(Work, BBBB_LENGTH_ABCI, IRREPX, 128+IOFF)
         Call Getall(Work, ABAB_LENGTH_ABCI, IRREPX, 129)
         Call Putall(Work, ABAB_LENGTH_ABCI, IRREPX, 129+IOFF)
         Call Getall(Work, ABBA_LENGTH_ABCI, IRREPX, 130)
         Call Putall(Work, ABBA_LENGTH_ABCI, IRREPX, 130+IOFF)

      Else
CSSS         AAAA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,127))
CSSS         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))

CSSS         Call Inipck(IRREPX,ISYTYP(1,127),ISYTYP(2,127),127+IOFF,IMODE,
CSSS     +               0,1)

         Call Inipck(IRREPX,ISYTYP(1,130),ISYTYP(2,130),130+IOFF,IMODE,
     +               0,1)

CSSS         Call Getall(Work, AAAA_LENGTH_ABCI, IRREPX, 127)
CSSS         Call Putall(Work, AAAA_LENGTH_ABCI, IRREPX, 127+IOFF)
         Call Getall(Work, ABBA_LENGTH_ABCI, IRREPX, 130)
         Call Putall(Work, ABBA_LENGTH_ABCI, IRREPX, 130+IOFF)

      Endif

C Hbar(IA,JK)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,107),ISYTYP(2,107))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,108),ISYTYP(2,108))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,109),ISYTYP(2,109))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))

         Call Inipck(IRREPX,ISYTYP(1,127),ISYTYP(2,107),107+IOFF,IMODE,
     +               0,1)
         Call Inipck(IRREPX,ISYTYP(1,128),ISYTYP(2,108),108+IOFF,IMODE,
     +               0,1)
         Call Inipck(IRREPX,ISYTYP(1,129),ISYTYP(2,109),109+IOFF,IMODE,
     +               0,1)
         Call Inipck(IRREPX,ISYTYP(1,130),ISYTYP(2,110),110+IOFF,IMODE,
     +               0,1)

         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 107)
         Call Putall(Work, AAAA_LENGTH_IJKA, IRREPX, 107+IOFF)
         Call Getall(Work, BBBB_LENGTH_IJKA, IRREPX, 108)
         Call Putall(Work, BBBB_LENGTH_IJKA, IRREPX, 108+IOFF)
         Call Getall(Work, ABAB_LENGTH_IJKA, IRREPX, 109)
         Call Putall(Work, ABAB_LENGTH_IJKA, IRREPX, 109+IOFF)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 110+IOFF)

      Else
CSSS         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,107),ISYTYP(2,107))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))

CSSS         Call Inipck(IRREPX,ISYTYP(1,107),ISYTYP(2,107),107+IOFF,IMODE,
CSSS     +               0,1)
         Call Inipck(IRREPX,ISYTYP(1,130),ISYTYP(2,110),110+IOFF,IMODE,
     +               0,1)

CSSS         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 107)
CSSS         Call Putall(Work, AAAA_LENGTH_IJKA, IRREPX, 107+IOFF)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)
         Call Putall(Work, ABBA_LENGTH_IJKA, IRREPX, 110+IOFF)

      Endif

      Return
      End


