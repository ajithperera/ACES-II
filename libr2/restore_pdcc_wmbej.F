      Subroutine restore_pdcc_wmbej(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Integer AAAA_LENGTH_MBEJ,BBBB_LENGTH_MBEJ,AABB_LENGTH_MBEJ
      Integer BBAA_LENGTH_MBEJ,ABAB_LENGTH_MBEJ,BABA_LENGTH_MBEJ
      Double Precision Mone, One

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
C             -
C The P(D) CC W(mb,ej) is needed for the L2INL2 ring contributions.
C These are computed in post_cc_mods are stored in 165-169 list.
C Restore them on 54-59 lists.
C
      IRREPX = 1

      AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
      AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
      ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

      If (Iuhf .NE. 0) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
      End if

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") "The checksums of Hbar(me,bj) at Entry:"
      Write(6,*)
      If (IUhf .Gt. 0) Then
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 164)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 165)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 166)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 167)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 168)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 169)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
      Else
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 164)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 166)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 168)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
      Endif

      If (IUhf .Gt. 0) Then
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
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, ABAB_LENGTH_MBEJ)
      Endif

#endif


      IF (IUhf .GT. 0) Then

         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 164)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 165)
         Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 166)
         Call putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 167)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 168)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 169)
         Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)

      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 164)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 166)
         Call putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 168)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)

      Endif

      Return
      End
 
