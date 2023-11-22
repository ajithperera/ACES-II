
      program main
      implicit none

c INTERNAL VARIABLES
c   o the number of doubles useable by this program
      integer dcoredim
c   o the MPI error code(s)
      integer iMPIErr

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "method.com"

c ----------------------------------------------------------------------

c   o initialize the MPI environment
      call MPI_INIT(iMPIErr)
      if (iMPIErr.ne.0) then
         print *, '@(XCC)MAIN: unable to initialize the MPI environment'
         call xcc_errex
      end if

c   o initialize the ACES environment but do not allocate memory
      call aces_init(icore,i0,0,iuhf,.false.)

c   o initialize xcc and allocate memory
      call xcc_init(dcoredim)

c   o mask the offset, the dimension, and the type
      call xcc_driver(icore(i0),dcoredim)

c   o finalize the ACES environment
      call aces_fin

c   o finalize the MPI environment
      call MPI_FINALIZE(iMPIErr)
      if (iMPIErr.ne.0) then
         print *, '@(XCC)MAIN: unable to finalize the MPI environment'
         call xcc_errex
      end if

c ----------------------------------------------------------------------

      call exit(0)
      end

