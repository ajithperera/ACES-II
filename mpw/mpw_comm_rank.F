
#include "mpwf.h"

      subroutine MPW_Comm_rank(comm,rank,ierror)
      implicit none

      integer comm, rank, ierror
      F_MPI_INT w_comm, w_rank, w_ierror

      w_comm = comm
      call MPI_Comm_rank(w_comm,w_rank,w_ierror)
      rank = w_rank
      ierror = w_ierror

      return
      end

