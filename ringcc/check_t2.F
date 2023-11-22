      Subroutine Check_t2(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"
#include "syminf.com"

      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)
  
#ifdef _NOSKIP
      Write(*,*) 
      Nmult = 1
      If (uhf) Nmult=2 
      do ispin = 1, Nmult 
      do irrep = 1, Nirrep
         nroww =Irpdpd(irrep,18+ispin)
         ncolw =Irpdpd(irrep,20+ispin)
         write(6,*) nroww,ncolw
         call getlst(work,1,ncolw,1,irrep,2+Ispin)
         Call Checksum("<AB||IJ>:",Work,nroww*ncolw) 
      enddo 
      enddo 
      Write(*,*) 
#endif 

C T2(IJ,AB)
C#ifdef _NOSKIP
C      do irrep = 1, Nirrep
C         nroww =Irpdpd(irrep,1)
C         ncolw =Irpdpd(irrep,3)
C         call getlst(work,1,ncolw,1,irrep,14)
C         Call Checksum("W:",Work,nroww*ncolw) 
C         nrowt =Irpdpd(irrep,19)
C         ncolt =Irpdpd(irrep,21)
C         call getlst(work,1,ncolt,1,irrep,44)
C         Call Checksum("T:",Work,nrowt*ncolt) 
C      enddo 

      Irrepx = 1
      Length_61=IDSYMSZ(IRREPX,ISYTYP(1,61),ISYTYP(2,61))
      Length_62=IDSYMSZ(IRREPX,ISYTYP(1,62),ISYTYP(2,62))
      Length_63=IDSYMSZ(IRREPX,ISYTYP(1,63),ISYTYP(2,63))

      Write(6,*)
      Call Getall(Work, Length_61, Irrepx, 61)
      Call checksum("T2-AAAA :",Work,Length_61)
      Call Getall(Work, Length_62, Irrepx, 62)
      Call checksum("T2-BBBB :",Work,Length_62)
      Call Getall(Work, Length_63, Irrepx, 63)
      Call checksum("T2-ABAB :",Work,Length_63)

      Write(*,*) 
      IRREPX = 1
      IF (UHF) Then
         AAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
         BBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45))
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         Call Getall(Work, AAAA_LENGTH_IJAB, IRREPX, 44)
         Call Checksum("T2AAAA_ABIJ", Work, AAAA_LENGTH_IJAB)
         Call Getall(Work, BBBB_LENGTH_IJAB, IRREPX, 45)
         Call Checksum("T2BBBB_abij", Work, BBBB_LENGTH_IJAB)
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 46)
         Call Checksum("T2ABAB_AbIj", Work, ABAB_LENGTH_IJAB)

         AAAA_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,34),ISYTYP(2,34))
         Call Getall(Work, AAAA_LENGTH_AJBI, IRREPX, 34)
         Call Checksum("T2AAAA_AJBI", Work, AAAA_LENGTH_AJBI)
         BBBB_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,35),ISYTYP(2,35))
         Call Getall(Work, BBBB_LENGTH_AJBI, IRREPX, 35)
         Call Checksum("T2BBBB_AJBI", Work, BBBB_LENGTH_AJBI)
         BBAA_LENGTH_BJAI = IDSYMSZ(IRREPX,ISYTYP(1,36),ISYTYP(2,36))
         Call Getall(Work, BBBB_LENGTH_BJAI, IRREPX, 36)
         Call Checksum("T2BBAA_BJAI", Work, BBAA_LENGTH_BJAI)
         AABB_LENGTH_AIBJ = IDSYMSZ(IRREPX,ISYTYP(1,37),ISYTYP(2,37))
         Call Getall(Work, BBBB_LENGTH_AIBJ, IRREPX, 37)
         Call Checksum("T2AABB_AIBJ", Work, AABB_LENGTH_AIBJ)
         BAAB_LENGTH_BIAJ = IDSYMSZ(IRREPX,ISYTYP(1,38),ISYTYP(2,38))
         Call Getall(Work, BAAB_LENGTH_BIAJ, IRREPX, 38)
         Call Checksum("T2BIAJ_BIAJ", Work, BAAB_LENGTH_BIAJ)
         ABBA_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,39),ISYTYP(2,39))
         Call Getall(Work, ABBA_LENGTH_AJBI, IRREPX, 39)
         Call Checksum("T2AJBI_ABBA", Work, ABBA_LENGTH_AJBI)

      Else
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 46)
         Call Checksum("T2ABAB_AbIj", Work, ABAB_LENGTH_IJAB)
         AAAA_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,34),ISYTYP(2,34))
         Call Getall(Work, AAAA_LENGTH_AJBI, IRREPX, 34)
         Call Checksum("T2AAAA_AJBI", Work, AAAA_LENGTH_AJBI)

         AABB_LENGTH_AIBJ = IDSYMSZ(IRREPX,ISYTYP(1,37),ISYTYP(2,37))
         Call Getall(Work, BBBB_LENGTH_AIBJ, IRREPX, 37)
         Call Checksum("T2AABB_AIBJ", Work, AABB_LENGTH_AIBJ)
         aBBA_LENGTH_AJBI = IDSYMSZ(IRREPX,ISYTYP(1,39),ISYTYP(2,39))
         Call Getall(Work, ABBA_LENGTH_AJBI, IRREPX, 39)
         Call Checksum("T2AJBI_ABBA", Work, ABBA_LENGTH_AJBI)
      Endif

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
C#endif   

      Write(6,*)
      Return
      End


