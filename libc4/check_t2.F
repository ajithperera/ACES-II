      Subroutine Check_T2(Work,Length,Iuhf,T2_fine)
      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Length)
      Double Precision Nm_aaaa,Nm_bbbb,Nm_abab,Nm_aa,Nm_bb
      Double Precision Nm_max
      Logical T2_fine

      Logical null 
      Data onem,One,Dzero /-1.0,1.0,0.0D0/
#include "sympop.com"
#include "syminf.com"

C T2-ABIJ
      Irrepx = 1
      Length_44=IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44))
      If (Iuhf .NE. 0) Length_45=IDSYMSZ(IRREPX,ISYTYP(1,45),
     +                                   ISYTYP(2,45))
      Length_46=IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46))

      Write(6,*) 
      Call Getall_c4(Work, Length_44, Irrepx, 44)
      Call checksum("T2-AAAA :",Work,Length_44,Nm_aaaa)
      If (Iuhf .ne. 0) Then
      Call Getall_c4(Work, Length_45, Irrepx, 45)
      Call checksum("T2-BBBB :",Work,Length_45,Nm_bbbb)
      Endif 
      Call Getall_c4(Work, Length_46, Irrepx, 46)
      Call checksum("T2-ABAB :",Work,Length_46,Nm_abab)

C T1-AI
      Length_aa =Irpdpd(Irrepx,9)
      call getlst_c4(Work,1,1,1,3,90)
      Call checksum("T1-AA   :",Work,Length_aa,Nm_aa)
      If (Iuhf .Ne. 0) Then 
      Length_bb =Irpdpd(Irrepx,10)
      call getlst_c4(Work,1,1,1,4,90)
      Call checksum("T1-BB   :",Work,Length_bb,Nm_bb)
      Endif 
      Write(6,*) 

      Nm_max = Max(Nm_aa,Nm_bb,Nm_aaaa,Nm_bbbb,Nm_abab)
      T2_fine = .False.
      If (NM_max .Ne. Dzero) T2_fine = .True.

      Return
      End

      subroutine checksum(sTag,a,n,nm)
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
      write(*,'(2a,3x,3e24.12)') '@CHECKSUM: ',sTag,sm,nm,dsqrt(nm)
      return
      end


