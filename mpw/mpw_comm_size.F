
#include "mpwf.h"

      subroutine MPW_Comm_size(comm,size,ierror)
      implicit none

      integer comm, size, ierror
      F_MPI_INT w_comm, w_size, w_ierror

      w_comm = comm
      call MPI_Comm_size(w_comm,w_size,w_ierror)
      size = w_size
      ierror = w_ierror

      return
      end

