
#include "mpwf.h"

      subroutine MPW_Initialized(flag,ierror)
      implicit none

      logical flag
      integer ierror
      F_MPI_LOG w_flag
      F_MPI_INT w_ierror

      call MPI_Initialized(w_flag,w_ierror)
      flag = w_flag
      ierror = w_ierror

      return
      end

