      Subroutine Checkgamma(Work,Maxcor,Iuhf)

      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF 
     
#include "machsp.com"
#include "sympop.com"
#include "parcc.com"

      UHF = .FALSE.
      UHF = (IUhf .EQ. 1)

      IRREPX = 1

C Gamma(AB,CD)

      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,111),ISYTYP(2,111))
         BBBB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,112),ISYTYP(2,112))
         ABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,113),ISYTYP(2,113))
         Call Getall(Work, AAAA_LENGTH_ABCD, IRREPX, 51)
         Call Checksum("AAAA_ABCD", Work, AAAA_LENGTH_ABCD)
         Call Getall(Work, BBBB_LENGTH_ABCD, IRREPX, 52)
         Call Checksum("BBBB_ABCD", Work, BBBB_LENGTH_ABCD)
         Call Getall(Work, ABAB_LENGTH_ABCD, IRREPX, 53)
         Call Checksum("ABAB_ABCD", Work, ABAB_LENGTH_ABCD)
      Else
         ABAB_LENGTH_ABCD = IDSYMSZ(IRREPX,ISYTYP(1,113),ISYTYP(2,113))
         Call Getall(Work, ABAB_LENGTH_ABCD, IRREPX, 113)
         Call Checksum("ABAB_ABCD", Work, ABAB_LENGTH_ABCD)
      Endif

C Gamma(MN,IJ)
  
      Write(*,*) 
      IF (UHF) Then
         AAAA_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,111),ISYTYP(2,111))
         BBBB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,112),ISYTYP(2,112))
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,113),ISYTYP(2,113))
         Call Getall(Work, AAAA_LENGTH_MNIJ, IRREPX, 51)
         Call Checksum("AAAA_MNIJ", Work, AAAA_LENGTH_MNIJ)
         Call Getall(Work, BBBB_LENGTH_MNIJ, IRREPX, 52)
         Call Checksum("BBBB_MNIJ", Work, BBBB_LENGTH_MNIJ)
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Checksum("ABAB_MNIJ", Work, ABAB_LENGTH_MNIJ)
      Else
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,113),ISYTYP(2,113))
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 113)
         Call Checksum("ABAB_MNIJ", Work, ABAB_LENGTH_MNIJ)
      Endif

C Gamma(AB,CI)

      Write(*,*) 
      IF (UHF) Then
         AAAA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,127))
         BBBB_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,128),ISYTYP(2,128))
         ABAB_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,129),ISYTYP(2,129))
         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))
         Call Getall(Work, AAAA_LENGTH_AIBC, IRREPX, 127)
         Call Checksum("AAAA_AIBC", Work, AAAA_LENGTH_AIBC)
         Call Getall(Work, BBBB_LENGTH_AIBC, IRREPX, 128)
         Call Checksum("BBBB_AIBC", Work, BBBB_LENGTH_AIBC)
         Call Getall(Work, ABAB_LENGTH_AIBC, IRREPX, 129)
         Call Checksum("ABAB_AIBC", Work, ABAB_LENGTH_AIBC)
         Call Getall(Work, ABBA_LENGTH_AIBC, IRREPX, 130) 
         Call Checksum("ABBA_AIBC", Work, ABBA_LENGTH_AIBC)
      Else
         ABBA_LENGTH_AIBC = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))
         Call Getall(Work, ABBA_LENGTH_AIBC, IRREPX, 130) 
         Call Checksum("ABBA_AIBC", Work, ABBA_LENGTH_AIBC)
      Endif 

C Gamma(IJ,KA)

      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,107),ISYTYP(2,107))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,108),ISYTYP(2,108))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,109),ISYTYP(2,109))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))
         Call Getall(Work, AAAA_LENGTH_IJKA, IRREPX, 107)
         Call Checksum("AAAA_IJKA", Work, AAAA_LENGTH_IJKA)
         Call Getall(Work, BBBB_LENGTH_IJKA, IRREPX, 108)
         Call Checksum("BBBB_IJKA", Work, BBBB_LENGTH_IJKA)
         Call Getall(Work, ABAB_LENGTH_AIBC, IRREPX, 109)
         Call Checksum("ABAB_IJKA", Work, ABAB_LENGTH_IJKA)
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)
         Call Checksum("ABBA_IJKA", Work, ABBA_LENGTH_IJKA)
      Else
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))
         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 110)
         Call Checksum("ABBA_IJKA", Work, ABBA_LENGTH_IJKA,S)
      Endif

C Gamma(MB,EJ)

      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,123),ISYTYP(2,123))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,124),ISYTYP(2,124))
         ABBA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,125),ISYTYP(2,125))
         BAAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,126),ISYTYP(2,126))
	 ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,118),ISYTYP(2,118))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,117),ISYTYP(2,117))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 123)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 124)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ)
         Call Getall(Work, ABBA_LENGTH_MBEJ, IRREPX, 125)
         Call Checksum("ABBA_MBEJ", Work, ABBA_LENGTH_MBEJ)
         Call Getall(Work, BAAB_LENGTH_MBEJ, IRREPX, 126)
         Call Checksum("BAAB_MBEJ", Work, BAAB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 118)
         Call Checksum("ABAB_MBEJ", Work, ABAB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 117)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,123),ISYTYP(2,123))
         ABBA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,125),ISYTYP(2,125))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,118),ISYTYP(2,118))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 123)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, ABBA_LENGTH_MBEJ, IRREPX, 125)
         Call Checksum("ABBA_MBEJ", Work, ABBA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 118)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
      Endif

C Gamma(IJ,AB)

      Write(*,*)
      IF (UHF) Then
         AAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,114),ISYTYP(2,114))
         BBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,115),ISYTYP(2,115))
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,116),ISYTYP(2,116))
         Call Getall(Work, AAAA_LENGTH_IJAB, IRREPX, 114)
         Call Checksum("AAAA_IJAB", Work, AAAA_LENGTH_IJAB)
         Call Getall(Work, BBBB_LENGTH_IJAB, IRREPX, 115)
         Call Checksum("BBBB_IJAB", Work, BBBB_LENGTH_IJAB)
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 116)
         Call Checksum("ABAB_IJAB", Work, ABAB_LENGTH_IJAB)
      Else
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,116),ISYTYP(2,116))
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 116)
         Call Checksum("ABAB_IJAB", Work, ABAB_LENGTH_IJAB)
      Endif

C DCC- Gamma(IJ,AB)
      If (Ispar .and. Coulomb) Then
      IF (UHF) Then
         AAAA_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,214),ISYTYP(2,214))
         BBBB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,215),ISYTYP(2,215))
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,216),ISYTYP(2,216))
         Call Getall(Work, AAAA_LENGTH_IJAB, IRREPX, 214)
         Call Checksum("AAAA_IJAB", Work, AAAA_LENGTH_IJAB)
         Call Getall(Work, BBBB_LENGTH_IJAB, IRREPX, 215)
         Call Checksum("BBBB_IJAB", Work, BBBB_LENGTH_IJAB)
         Call Getall(Work, ABAB_LENGTH_MNIJ, IRREPX, 216)
         Call Checksum("ABAB_IJAB", Work, ABAB_LENGTH_IJAB)
      Else
         ABAB_LENGTH_IJAB = IDSYMSZ(IRREPX,ISYTYP(1,216),ISYTYP(2,216))
         Call Getall(Work, ABAB_LENGTH_IJAB, IRREPX, 216)
         Call Checksum("ABAB_IJAB", Work, ABAB_LENGTH_IJAB)
      Endif
      Endif 
      Write(6,*)
      Return
      End

c This routine prints the sum of elements a(1:n) and their L2 vector norm.
c If two vectors produce the same sum but a different norm, then they
c have the same contributions in different distributions. If element
c order is to be checked, then use checksum2.

      subroutine checksum(sTag,a,n)
      implicit none

      character*(*) sTag
      integer n
      double precision a(n)

      double precision sm, nm, s
      integer i

      sm=0.d0
      nm=0.d0
      do i=1,n
         sm=sm+a(i)
         nm=nm+a(i)*a(i)
      end do
      write(*,'(2a,3x,2e24.12)') '@CHECKSUM: ',sTag,sm,nm

      s= sqrt(nm)
      
      return
      end

