
#include "mpwf.h"

      subroutine MPW_Allreduce(sendbuf,recvbuf,
     &                         count,datatype,op,comm,ierror)
      implicit none

      integer sendbuf, recvbuf
      integer count, datatype, op, comm, ierror
      F_MPI_INT w_count, w_datatype, w_op, w_comm, w_ierror
#include "mpif.h"

      w_count = count
      w_datatype = datatype
      w_op = op
      w_comm = comm

#ifdef F_64BIT /* 64-bit ACES */
#   ifndef _64BIT_MPI /* 32-bit MPI */
      if (datatype.eq.MPI_INTEGER) w_datatype = MPI_INTEGER8
#   endif
#else /* 32-bit ACES */
#   ifdef _64BIT_MPI /* 64-bit MPI */
      if (datatype.eq.MPI_INTEGER) w_datatype = MPI_INTEGER4
#   endif
#endif

      call MPI_Allreduce(sendbuf,recvbuf,
     &                   w_count,w_datatype,w_op,w_comm,w_ierror)
      ierror = w_ierror

      return
      end

