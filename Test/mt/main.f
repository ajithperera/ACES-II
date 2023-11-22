
c#define _DATA










      program main
      implicit none

      double precision alpha, a(65536,5), b(5,500),
     &                 beta,  c(65536,500)
      double precision cksm
      integer i, j

      data alpha, beta /1.d0, 0.d0/

      do j = 1, 5
      do i = 1, 65536
         a(i,j) = 1.d0
      end do
      end do
      do j = 1, 500
      do i = 1, 5
         b(i,j) = 2.d0
      end do
      end do
      do j = 1, 500
      do i = 1, 65536
         c(i,j) = 0.d0
      end do
      end do


      call dgemm('N','N',65536,500,5,
     &           alpha, a, 65536,
     &                  b, 5,
     &           beta,  c, 65536)

      cksm = 0.d0
      do j = 1, 500
      do i = 1, 65536
         cksm = cksm + c(i,j)
      end do
      end do
      print *, cksm

      end

