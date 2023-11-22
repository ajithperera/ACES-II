      Subroutine Check_hbar(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"
#include "syminf.com"

      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)
     
      Irrepx = 1
      Write(6,*)
      IF (UHF) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
      Endif

      Write(6,*)
      Return
      End


