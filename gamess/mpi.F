c----------- MPI wrapper routines ------------------------------------
c---------------------------------------------------------------------
      subroutine MPIAR1i (ni,n)
c computes n = \sum_{i=1}^{np} n_i
      integer ni,n
      integer iErr
#include "mpif.h"
#ifdef _HAVE_MPI
      call MPW_Allreduce(ni,n,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,iErr)
      if (iErr.ne.MPI_SUCCESS) then
         print *, '@MPIAR1i: Error from MPI_Allreduce'
         call aces_exit(1)
      end if
#else
      n = ni
#endif
      return
      end
c-------------------------------------------------------------
      subroutine MPIAR (xmat1,xmat2,n)
      implicit double precision (a-h,o-z)
      dimension xmat1(*), xmat2(*)
      parameter (maxl=4000000)  ! max buffer length possible in reduce
c                               ! On or O2K it empirically turned out
c                               ! to be something like 4,000,000
#include "mpif.h"

#ifdef _HAVE_MPI
c the lines below should not be needed if the work array XMAT2 were always
c available;  however, we don't want to allocate memory (in the main module)
c for XINT2 in sequential runs; since we restrict the use of USEMPI to 
c the MPI wrapper routines, the main module defines sequential runs as those
c with nprocs=1.
      call getnp (nprocs)
      if (nprocs.eq.1) return

c Find nn - the number of full-length batches
      nn = n/maxl
c Reduce the full batches
      do i=1,nn
       ii = (i-1)*maxl + 1
       call MPW_Allreduce(xmat1(ii),xmat2(ii),maxl,MPI_DOUBLE_PRECISION,
     >                    MPI_SUM,MPI_COMM_WORLD,iErr)
      if(iErr.ne.MPI_SUCCESS) then
       write(*,*)'Error in MPIAR in loop',i,ii,maxl
       call aces_exit(1)
      endif
      end do
c Reduce the last batch:
      ii = nn*maxl+1
      nnn = n - nn*maxl
      call MPW_Allreduce(xmat1(ii),xmat2(ii),nnn,MPI_DOUBLE_PRECISION,
     >                   MPI_SUM,MPI_COMM_WORLD,iErr)
      if(iErr.ne.MPI_SUCCESS) then
       write(*,*)'Error in MPIAR in last batch',ii,nnn
       call aces_exit(1)
      endif
c
      call XCOPY (n,xmat2,1,xmat1,1)
#endif

      end
