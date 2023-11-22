
      program main
      implicit none

c EXTERNAL FUNCTIONS
      integer cd_nodename

c INTERNAL VARIABLES
      integer iErr

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "method.com"
#include "machsp.com" /* for iintfp */

c ----------------------------------------------------------------------

c   o initialize the MPI environment
      call MPI_Init(iErr)
      if (iErr.ne.0) then
         print *, '@MAIN: unable to initialize the MPI environment'
         call aces_exit(iErr)
      end if

c   o change to the temp directory
      iErr = cd_nodename()
      if (iErr.ne.0) then
         print *, '@MAIN: unable to change to the temp directory'
         print *, '       attempting to run in current directory'
      end if

c   o initialize the ACES environment and allocate memory
      call aces_init(icore,i0,icrsiz,iuhf,.true.)

c   o mask the offset, the dimension, and the type
      call ft2w(icore(i0),icrsiz/iintfp)

c   o finalize the ACES environment
      call aces_fin

c   o finalize the MPI environment
      call MPI_Finalize(iErr)
      if (iErr.ne.0) then
         print *, '@MAIN: unable to finalize the MPI environment'
         call aces_exit(iErr)
      end if

c ----------------------------------------------------------------------

      call exit(0)
c     end program main
      end

