      Subroutine restore_cc_wmbej(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Integer AAAA_LENGTH_MBEJ,BBBB_LENGTH_MBEJ,AABB_LENGTH_MBEJ
      Integer BBAA_LENGTH_MBEJ,ABAB_LENGTH_MBEJ,BABA_LENGTH_MBEJ
      Double Precision Mone, One

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
C      
C CC W(mb,ej) is needed for L1 contribution from W(mb,ej). This 
C is computed in Post_vcc_mods and stored in 256-269. Restore
C them to 54-59 before the L1inl1 contractions. 
C
      IRREPX = 1

      IF (IUhf .GT. 0) Then

         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 264)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 265)
         Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 266)
         Call putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 267)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 268)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 269)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 59)

      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 264)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 266)
         Call putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 268)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)

      Endif

      Return
      End
 
