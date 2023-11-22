
#include "mpwf.h"

      subroutine MPW_Finalize(ierror)
      implicit none

      integer ierror
      F_MPI_INT w_ierror

      call MPI_Finalize(w_ierror)
      ierror = w_ierror

      return
      end

