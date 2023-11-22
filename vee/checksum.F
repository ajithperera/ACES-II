c This routine prints the sum of elements a(1:n) and their L2 vector norm.
c If two vectors produce the same sum but a different norm, then they
c have the same contributions in different distributions. If element
c order is to be checked, then use checksum2.

      subroutine checksum(sTag,a,n,s)
      implicit none

      character*(*) sTag
      integer n
      double precision a(n)

      double precision sm, nm, s
      integer i

CSSS      return 
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

