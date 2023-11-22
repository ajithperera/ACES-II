
#include "mpwf.h"

      subroutine MPW_Get_processor_name(name,resultlen,ierror)
      implicit none

      character*(*) name
      integer resultlen, ierror
      F_MPI_INT w_resultlen, w_ierror

      call MPI_Get_processor_name(name,w_resultlen,w_ierror)
      resultlen = w_resultlen
      ierror = w_ierror

      return
      end

