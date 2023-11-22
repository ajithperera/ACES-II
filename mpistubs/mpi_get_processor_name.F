
      subroutine mpi_get_processor_name(name,resultlen,ierror)
      character*(*) name
      integer resultlen, ierror
#include "mpif.h"
      name = 'localhost'
      resultlen = 9
      ierror = MPI_SUCCESS
      return
      end

