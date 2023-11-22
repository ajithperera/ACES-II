
c This routine drives the calculation.

c INPUT
c    integer n : the number of doubles available at dcore(1)
c    double precision dcore(n) : the core space for doubles

      subroutine xcc_driver(n,dcore)
      implicit none

c ARGUMENT LIST
      integer n
      double precision dcore(n)

c INTERNAL VARIABLES
      integer i, ierr

c EXTERNAL FUNCTIONS
      integer get_t

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c   o initialize dcore, the pointer to free memory, and the overlap
      dcoresize=n
      do i = 1,dcoresize
         dcore(i)=0.0D0
      end do
      freecore=1
      overlap=1.0D0

c   o load the T1 and T2 amplitudes into dcore and compute their dot
c     products
      ierr=get_t(dcore)
      if (ierr.ne.0) then
         write(*,*) '@XCC: The T amplitudes could not be loaded.'
         call errex
         stop 1
      end if

c   o calculate the 3rd order overlap term(s) <0|adj(T2)*T1^2|0>
      call ovrlp_t1t1t2(dcore)

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

