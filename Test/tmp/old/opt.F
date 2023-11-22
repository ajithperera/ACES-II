
#include "aces.h"
#define LOOP_MAX 10000
#define FUNC 500
#define DIM 50000
#define COL 1000
      program main
      implicit none
      double precision src(DIM+1), dest(DIM+1), ddot, dtmp, dsum
      integer i, j, k, len_row, len_col, iadd, ioff
      parameter (len_row=COL, len_col=DIM/COL)
#include "timeinfo.com"
      call timer(0)

      do i = 1, DIM
         src(i)  = i
         dest(i) = 1/i
      end do

      dsum = 0.
      dtmp = -1.
      call timer(1)
      do k = 1, LOOP_MAX
         ddot = 0.
         do i = 1, DIM
            ddot = ddot + ( src(i) * dest(i) )
         end do
         dsum = dsum + ( dtmp * ddot )
      end do
      call timer(1)
      print *, timenew, dsum

      do i = 1, DIM
         src(i)  = i
         dest(i) = 1/i
      end do

      dsum = 0.
      dtmp = 1.
      call timer(1)
      do k = 1, LOOP_MAX
         ddot = 0.
         do i = 1, DIM
            ddot = ddot + ( src(i) * dest(i) )
         end do
         dsum = dsum + ( dtmp * ddot )
      end do
      dsum = -dsum
      call timer(1)
      print *, timenew, dsum

#ifdef _FAIL_COMPILE
      k = 0
      call timer(1)
      do i = 1, LOOP_MAX
         j = FUNC
         do_j = .false.
         if (j.lt.i) do_j = (0.ne.j)
         if (j.eq.i) do_j = (1.ne.j)
         if (do_j) k = k + 1
      end do
      call timer(1)
      print *, timenew, k

      k = 0
      call timer(1)
      do i = 1, LOOP_MAX
         j = FUNC
         if (j.lt.i) then
            do_j = (0.ne.j)
         else
            if (j.eq.i) then
               do_j = (1.ne.j)
            else
               do_j = .false.
            end if
         end if
         if (do_j) k = k + 1
      end do
      call timer(1)
      print *, timenew, k
#endif /* _FAIL_COMPILE */

      end

#ifdef _FAIL_COMPILE

#define LEN 2
#define MAX_P 2
#define MAX_Q 3
#define MAX_R 4
#define ALPHA 1.0d0
#define FROM_Q 1
#define TO_Q 2

      program main
      implicit none

      integer    len,     p,       q,       r
      parameter (len=LEN, p=MAX_P, q=MAX_Q, r=MAX_R)
      integer    from_q,        to_q
      parameter (from_q=FROM_Q, to_q=TO_Q)
      double precision alpha
      parameter (alpha=ALPHA)

#define OMM_XX_XX_DDOT omm_qp_rq_ddot

      integer    a_row,   a_col
      parameter (a_row=q, a_col=p)
      double precision gea(len,a_row,a_col)

      integer    b_row,   b_col
      parameter (b_row=p, b_col=r)
      double precision geb(len,b_row,b_col)

      integer    c_row,   c_col
      parameter (c_row=p, c_col=r)
      double precision c(c_row,c_col)

      double precision sya(len,c_col*(c_col-1)/2)

      integer i, j, l

      print *, 'len = ', len
      print *, 'p = ', p
      print *, 'q = ', q
      print *, 'r = ', r
c      print *, 'from_q = ', from_q
c      print *, 'to_q = ', to_q
      print *, 'alpha = ', alpha

c   o init c(*)
      do j = 1, c_col
      do i = 1, c_row
         c(i,j) = 0.0d0
      end do
      end do

c   o init gea(*)
      do j = 1, a_col
      do i = 1, a_row
      do l = 1, len
         gea(l,i,j) = l*(i)
      end do
c         print *, 'A(*,',i,',',j,'): ',(gea(l,i,j),l=1,len)
      end do
      end do

c   o init geb(*)
      do j = 1, b_col
      do i = 1, b_row
      do l = 1, len
         geb(l,i,j) = l*(j)
      end do
         print *, 'B(*,',i,',',j,'): ',(geb(l,i,j),l=1,len)
      end do
      end do

c   o init sya(*)
      do j = 1, c_col-1
      do i = 0, j-1
      do l = 1, len
         sya(l,1+i+j*(j-1)/2) = l*(i)
      end do
         print *, 'A(*,',1+i+j*(j-1)/2,'): ',
     &            (sya(l,1+i+j*(j-1)/2),l=1,len)
      end do
      end do

      call omm_pq_qltq_ddot(len,b_row,b_col,
     &                      alpha,geb,b_row,0,b_col-1,
     &                            sya,
     &                            c,c_row)
c      call OMM_XX_XX_DDOT(len,p,r,q,
c     &                    alpha,gea,a_row,
c     &                          geb,b_row,
c     &                          c,  c_row)

      print *, 'C:'
      do i = 1, c_row
         print *, (c(i,j),j=1,c_col)
      end do

      end

#endif /* _FAIL_COMPILE */

