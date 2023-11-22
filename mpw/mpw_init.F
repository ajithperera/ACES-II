
#include "mpwf.h"

      subroutine MPW_Init(ierror)
      implicit none

      integer ierror
      F_MPI_INT w_ierror

      call MPI_Init(w_ierror)
      ierror = w_ierror

      return
      end

