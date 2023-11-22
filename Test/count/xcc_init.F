
c This routine:
c  o initializes the ACES environment (crapsi)
c  o allocates memory based on a suggestion
c  o initializes the silly I/O cache (inherited from crapsi)
c  o initializes the list information (inherited from crapsi)

c OUTPUT
c    integer dcoredim : the number of useable doubles from icore(i0)
c              BEWARE : Although this value is technically the same as
c                       'dcoresize' in the xcc common block, dcoredim
c                       is returned in order to partition the variables
c                       used by xcc_driver, xcc_init, and their
c                       respective dependencies.

      subroutine xcc_init(dcoredim)
      implicit none

c ARGUMENT LIST
      integer dcoredim

c INTERNAL VARIABLES
c   o the total suggested size for the icore array
      integer icoredim

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "machsp.com"
#include "method.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c   o initialize the ACES environment but do not allocate memory
      call crapsi(icore,iuhf,-1)

c   o size the problem
      call xcc_measure(icoredim,dcoredim)

c   o allocate the memory
      call aces_getmem(icoredim,icore,i0)
      if (i0.eq.-1) then
         write(*,*)
     &      '@XCC: Unable to allocate ',rshift(icoredim*iintln,20),
     &      ' MB needed for in-core amplitude processing.'
         call errex
         stop 1
      end if

c   o initialize the silly file cache
      call popcorn(icoredim)

c   o stupid crapsi . . . we have to do all this by hand if ientry=-1
      call getstf(icore(i0),iuhf,1)

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

