
c This routine:
c  o allocates memory based on a suggestion from xcc_measure()
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
c   o the total suggested size for the core array (in integers)
      integer icoredim

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "machsp.com"
#include "flags.com" /* for iflags(44) */
#include "method.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c   o process any command-line arguments
      call xcc_cl_proc

c   o size the problem
      call xcc_measure(icoredim,dcoredim)
      if (icoredim.lt.1) then
         print *, '@XCC_INIT: Failed to size memory.'
         call xcc_errex
      end if

c   o allocate the memory
#ifdef _DEBUG_XCC
      print *, '@XCC_INIT: Attempting to allocate ',
     &         icoredim,' integer-words...\n'
#endif /* _DEBUG_XCC */
      call aces_malloc(icoredim,icore,i0)
      if (icore(1).eq.0) then
         print *, '@XCC_INIT: Unable to allocate ',
     &            rshift(icoredim*iintln,20),
     &            ' MB needed for in-core amplitude processing.'
         call xcc_errex
      end if

c   o initialize the ACES I/O environment
      call aces_io_init(icore,i0,icoredim,.true.)

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

