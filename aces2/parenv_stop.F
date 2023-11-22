
      subroutine parenv_stop
      implicit none

      integer iErr

#include "parallel_aces.com"

#ifdef _HAVE_MPI
c   o finalize the MPI environment
      call MPW_Finalize(iErr)
      if (iErr.ne.0) then
         print *, '@MAIN: unable to finalize the MPI environment'
         call aces_exit(1)
      end if
#endif

      return
      end

