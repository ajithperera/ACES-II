
c This routine is the C1 equivalent of ovrlp_t1t1t2aa.F. It adds the
c quantity <0|adj(T2AA)*T1A^2|0> and its complex conjugate to the overlap.
c It can also be used for the T2BB piece if the argument ispin is 2;
c otherwise, if iuhf is 0, it doubles the T2AA term.

c INPUT
c    double precision dcore(*) : core space for doubles
c    integer          ispin    : Alpha/Beta signal (1=Alpha,2=Beta)

#include "disttype.h"

      subroutine nosym(dcore,ispin)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)
      integer ispin

c INTERNAL VARIABLES
      double precision dtmp, ovrlp_tmp
      integer t2_start, t1_start, ioff, nocc, nvrt
      integer i, j, a, b

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      nocc = pop(1,ispin)
      if (nocc.lt.2) return
      nvrt = vrt(1,ispin)
      if (nvrt.lt.2) return

      ioff = 0
      t2_start = ndx_t2(1,ispin)

      do j = 1,nocc-1
      do i = 0,j-1
      do b = 1,nvrt-1
      do a = 0,b-1
         write(*,'(1F17.12,3x,4I5)')
     &      dcore(t2_start+ioff),
     &      a+1, b+1, i+1, j+1
         ioff = ioff + 1
      end do
      end do
      end do
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

