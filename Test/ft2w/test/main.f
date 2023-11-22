
      program main
      implicit none

#ifdef _AIX
      include "mpif.h"
#else
      include "/usr/local/include/mpif.h"
#endif

      integer max_cpus
      parameter (max_cpus=4)

      integer irank, nprocs
      integer recvcounts(max_cpus), displs(max_cpus)
      integer sbuf(5), rbuf(5,max_cpus)

      integer i, j, iErr

c ----------------------------------------------------------------------

      call MPI_Init(iErr)
      if (iErr.ne.0) call exit(iErr)
      call MPI_Comm_size(MPI_COMM_WORLD,nprocs,iErr)
      if (iErr.ne.0) call exit(iErr)
      call MPI_Comm_rank(MPI_COMM_WORLD,irank,iErr)
      if (iErr.ne.0) call exit(iErr)

      do i = 1, max_cpus
         recvcounts(i) = 5
         displs(i)     = (i-1)*5
      end do
      do i = 1, 5
         rbuf(i,1+irank) = irank+1
      end do

      call MPI_Allgatherv(rbuf(1,1+irank), 5, MPI_INTEGER,
     &                    rbuf, recvcounts, displs, MPI_INTEGER,
     &                    MPI_COMM_WORLD, iErr)
      if (iErr.ne.0) call exit(iErr)
cint MPI_Allgatherv ( sendbuf, sendcount,  sendtype, 
c                     recvbuf, recvcounts, displs,   recvtype, comm )
cvoid             *sendbuf;
cint               sendcount;
cMPI_Datatype      sendtype;
cvoid             *recvbuf;
cint              *recvcounts;
cint              *displs;
cMPI_Datatype      recvtype;
cMPI_Comm          comm;

      do i = 1, 5
         print *, irank,':',(rbuf(i,j),j=1,nprocs)
      end do
      print '(/)'

      call MPI_Finalize(iErr)
      if (iErr.ne.0) call exit(iErr)

      end

