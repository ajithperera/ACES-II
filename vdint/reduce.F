
c This routine destroys column I of an array A(N,M) by shifting back the
c column range I+1 through M.

      subroutine reduce(a,n,i,m)
      implicit none
      integer n, i, m
      double precision a(n,m)
#include "machsp.com" /* for ifltln */
#ifdef _ASSERT
      if ((n.lt.1).or.(m.lt.i).or.(i.lt.1)) then
         print *, '@REDUCE: Assertion failed.'
         print *, '   rm col index = ',i
         print *, '   rows = ',n
         print *, '   cols = ',m
         call aces_exit(1)
      end if
#endif
      if (i.eq.m) return
      call c_memmove(a(1,i),a(1,i+1),n*(m-i)*ifltln)
      return
      end

