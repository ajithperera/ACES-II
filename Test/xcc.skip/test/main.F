
#define _DIM 100

      program main
      implicit none

      double precision darr(_DIM), dlsum, dsum
      integer          iarr(_DIM), ilsum, isum
      integer ntasks, myid, iErr
      integer i
#include </usr/local/include/mpif.h>

      call MPI_INIT(iErr)
      if (iErr.ne.MPI_SUCCESS) call mpi_errex('MPI_INIT',iErr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,iErr)
      if (iErr.ne.MPI_SUCCESS) call mpi_errex('MPI_COMM_SIZE',iErr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,iErr)
      if (iErr.ne.MPI_SUCCESS) call mpi_errex('MPI_COMM_RANK',iErr)

      if (myid.eq.0) print *, 'ntasks = ',ntasks
      dlsum = 0.d0
      ilsum = 0
      do i = (_DIM*myid)+1, _DIM*(myid+1)
         ilsum = ilsum + i
         dlsum = dlsum + i*1.d0
      end do
      print *, myid,': ilsum = ',ilsum,', dlsum = ',dlsum

      call MPI_REDUCE(ilsum,isum,1,
     &                MPI_INTEGER,MPI_SUM,
     &                0,MPI_COMM_WORLD,iErr)
      if (iErr.ne.MPI_SUCCESS) call mpi_errex('MPI_REDUCE',iErr)
      call MPI_REDUCE(dlsum,dsum,1,
     &                MPI_DOUBLE_PRECISION,MPI_SUM,
     &                0,MPI_COMM_WORLD,iErr)
      if (iErr.ne.MPI_SUCCESS) call mpi_errex('MPI_REDUCE',iErr)
      if (myid.eq.0) print *, '0: isum = ',isum,', dsum = ',dsum

      call MPI_FINALIZE(iErr)
      end

      subroutine mpi_errex(szFunc,iErr)
      implicit none
      character*(*) szFunc
      integer iErr
      print *, '@',szFunc,': error code = ',iErr
      call exit(1)
      end

