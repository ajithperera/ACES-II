
c#define _DATA

#define M	65536
#define LDA	65536
#define LDC	65536

#define K	5
#define LDB	5

#define N	500

      program main
      implicit none

      double precision alpha, a(LDA,K), b(LDB,N),
     &                 beta,  c(LDC,N)
      double precision cksm
      integer i, j

      data alpha, beta /1.d0, 0.d0/
#ifdef _DATA
      integer na, nb
      parameter (na=LDA*K, nb=LDB*N)
      data a, b /na*1.d0, nb*2.d0/
#else
      do j = 1, K
      do i = 1, LDA
         a(i,j) = 1.d0
      end do
      end do
      do j = 1, N
      do i = 1, LDB
         b(i,j) = 2.d0
      end do
      end do
      do j = 1, N
      do i = 1, LDC
         c(i,j) = 0.d0
      end do
      end do
#endif

      call dgemm('N','N',M,N,K,
     &           alpha, a, LDA,
     &                  b, LDB,
     &           beta,  c, LDC)

      cksm = 0.d0
      do j = 1, N
      do i = 1, LDC
         cksm = cksm + c(i,j)
      end do
      end do
      print *, cksm

      end

