
c This routine is the common error handler for MPI routines.

      subroutine mpi_errex(szFunc,iErr)
      implicit none

c ARGUMENTS
      character*(*) szFunc
      integer iErr

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      print *, '@',szFunc,': error code = ',iErr
      call xcc_errex

      end

