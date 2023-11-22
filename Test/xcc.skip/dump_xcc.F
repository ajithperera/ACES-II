
c This routine prints the contents of the xcc.com common block.

      subroutine dump_xcc
      implicit none

c INTERNAL VARIABLES
      integer i, j

c COMMON BLOCKS
#include "xcc.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      print *, '@DUMP_XCC: THE CONTENTS OF THE XCC COMMON BLOCK:'

      print *, '   overlap = ',overlap
      print *, '   <0|adj(T1)*T1|0> = ',t_sqr(1)
      print *, '   <0|adj(T2)*T2|0> = ',t_sqr(2)
      print *, '   <0|adj(T3)*T3|0> = ',t_sqr(3)

      print *, '   freecore = ',freecore,' of ',dcoresize

 2    format(a,i1,a,4i10)
      do i = 1, nirrep
         print 2, '   T1(',i,') @ ',(ndx_t1(i,j),j=1,2)
      end do
      do i = 1, nirrep
         print 2, '   T2(',i,') @ ',(ndx_t2(i,j),j=1,3)
      end do
      do i = 1, nirrep
         print 2, '   T3(',i,') @ ',(ndx_t3(i,j),j=1,4)
      end do
      do i = 1, nirrep
         print 2, ' T2p1(',i,') @ ',(ndx_t2p1(i,j),j=1,3)
      end do
      do i = 1, nirrep
         print 2, '   S1(',i,') @ ',(ndx_s1vo(i,j),j=1,2)
      end do
      print '(a,4a10)', '           ','I','A','i','a'
      do i = 1, nirrep
         print 2, '  eps(',i,') @ ',
     &                               ndx_eval(i,1,1),
     &                               ndx_eval(i,2,1),
     &                               ndx_eval(i,1,2),
     &                               ndx_eval(i,2,2)
      end do

c PRINT OUT off_pltp AND off_pq

      print *, '@DUMP_XCC: -------------------------------------'

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

