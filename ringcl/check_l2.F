      Subroutine Check_l2(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"
#include "syminf.com"

      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)

C T2(IJ,AB)
  
      Write(*,*) 
      IRREPX = 1
      IF (UHF) Then
         AAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,144),ISYTYP(2,144))
         BBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,145),ISYTYP(2,145))
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,146),ISYTYP(2,146))
         Call Getall(Work, AAAA_LENGTH_IJAB, IRREPX, 144)
         Call Checksum("L2AAAA_ABIJ:144", Work, AAAA_LENGTH_IJAB)
         Call Getall(Work, BBBB_LENGTH_IJAB, IRREPX, 145)
         Call Checksum("L2BBBB_abij:145", Work, BBBB_LENGTH_IJAB)
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 146)
         Call Checksum("L2ABAB_AbIj:146", Work, ABAB_LENGTH_IJAB)

         Call Getall(Work, AAAA_LENGTH_IJAB, IRREPX, 61)
         Call Checksum("L2AAAA_ABIJ:61", Work, AAAA_LENGTH_IJAB)
         Call Getall(Work, BBBB_LENGTH_IJAB, IRREPX, 62)
         Call Checksum("L2BBBB_abij:62", Work, BBBB_LENGTH_IJAB)
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 63)
         Call Checksum("L2ABAB_AbIj:63", Work, ABAB_LENGTH_IJAB)

         AAAA_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,134),ISYTYP(2,134))
         Call Getall(Work, AAAA_LENGTH_AJBI, IRREPX, 134)
         Call Checksum("L2AAAA_AJBI", Work, AAAA_LENGTH_AJBI)
         BBBB_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,135),ISYTYP(2,135))
         Call Getall(Work, BBBB_LENGTH_AJBI, IRREPX, 135)
         Call Checksum("L2BBBB_AJBI", Work, BBBB_LENGTH_AJBI)
         BBAA_LENGTH_BJAI = IDSYMSZ(IRREPX,ISYTYP(1,136),ISYTYP(2,136))
         Call Getall(Work, BBBB_LENGTH_BJAI, IRREPX, 136)
         Call Checksum("L2BBAA_BJAI", Work, BBAA_LENGTH_BJAI)
         AABB_LENGTH_AIBJ = IDSYMSZ(IRREPX,ISYTYP(1,137),ISYTYP(2,137))
         Call Getall(Work, BBBB_LENGTH_AIBJ, IRREPX, 137)
         Call Checksum("L2AABB_AIBJ", Work, AABB_LENGTH_AIBJ)
         BAAB_LENGTH_BIAJ = IDSYMSZ(IRREPX,ISYTYP(1,138),ISYTYP(2,138))
         Call Getall(Work, BAAB_LENGTH_BIAJ, IRREPX, 138)
         Call Checksum("L2BIAJ_BIAJ", Work, BAAB_LENGTH_BIAJ)
         ABBA_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,139),ISYTYP(2,139))
         Call Getall(Work, ABBA_LENGTH_AJBI, IRREPX, 139)
         Call Checksum("L2AJBI_ABBA", Work, ABBA_LENGTH_AJBI)

      Else
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,146),ISYTYP(2,146))
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 146)
         Call Checksum("L2ABAB_AbIj:146", Work, ABAB_LENGTH_IJAB)
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 63)
         Call Checksum("L2ABAB_AbIj:63", Work, ABAB_LENGTH_IJAB)

C         AAAA_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,134),ISYTYP(2,134))
C         Call Getall(Work, AAAA_LENGTH_AJBI, IRREPX, 134)
C         Call Checksum("L2AAAA_AJBI", Work, AAAA_LENGTH_AJBI)
C         AABB_LENGTH_AIBJ = IDSYMSZ(IRREPX,ISYTYP(1,137),ISYTYP(2,137))
C         Call Getall(Work, BBBB_LENGTH_AIBJ, IRREPX, 137)
C         Call Checksum("L2AABB_AIBJ", Work, AABB_LENGTH_AIBJ)
C         aBBA_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,139),ISYTYP(2,139))
C         Call Getall(Work, ABBA_LENGTH_AJBI, IRREPX, 139)
C         Call Checksum("L2AJBI_ABBA", Work, ABBA_LENGTH_AJBI)
      Endif
      Return

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


