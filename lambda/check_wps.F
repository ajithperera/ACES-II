        Subroutine Check_wps(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Integer AAAA_LENGTH_MBEJ,BBBB_LENGTH_MBEJ,AABB_LENGTH_MBEJ
      Integer BBAA_LENGTH_MBEJ,ABAB_LENGTH_MBEJ,BABA_LENGTH_MBEJ
      Integer AA_FAE_LENGTH, BB_FAE_LENGTH

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"

      Irrepx = 1
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*)     "                _         "
      Write(6,"(a)") "The checksums of W for for DCC:"
      Write(6,*)
      IF (Iuhf .Gt. 0) Then
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
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX,59)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
         Write(6,*)
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
         Write(6,*)
      Endif

      Write(6,*)
      Write(6,*)     "                =         "
      Write(6,"(a)") "The checksums of W for for DCC:"
      Write(6,*)
      IF (Iuhf .Gt. 0) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,123),ISYTYP(2,123))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,124),ISYTYP(2,124))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,118),ISYTYP(2,118))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,117),ISYTYP(2,117))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,125),ISYTYP(2,125))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,126),ISYTYP(2,126))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 123)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 124)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 118)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 117)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 125)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 126)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
         Write(6,*)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,123),ISYTYP(2,123))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,118),ISYTYP(2,118))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,125),ISYTYP(2,125))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 123)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 118)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 125)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Write(6,*)
      Endif

#endif
      Return
      End

