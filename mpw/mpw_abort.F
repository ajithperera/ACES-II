
#include "mpwf.h"

      subroutine MPW_Abort(comm,errorcode,ierror)
      implicit none

      integer comm, errorcode, ierror
      F_MPI_INT w_comm, w_errorcode, w_ierror

      w_comm = comm
      w_errorcode = errorcode
      call MPI_Abort(w_comm,w_errorcode,w_ierror)
      ierror = w_ierror

      return
      end

