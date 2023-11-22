
      subroutine mpi_initialized(flag,ierror)
      logical flag
      integer ierror
#include "mpif.h"
      flag = .true.
      ierror = MPI_SUCCESS
      return
      end

