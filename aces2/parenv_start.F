
      subroutine parenv_start
      implicit none

      integer iErr

#include "aces2.h"
#include "parallel_aces.com"

#ifdef _HAVE_MPI

c   o initialize the MPI environment
      call MPW_Init(iErr)
      if (iErr.ne.0) then
         print *,
     &         '@PARENV_START: unable to initialize the MPI environment'
         call aces_exit(1)
      end if

c   o initialize the parallel_aces common block
      call aces_com_parallel_aces

c   o change to the temp directory
      call cd_noderank(irank,iErr)
      if (iErr.ne.0) then
         print *, '@PARENV_START: unable to cd to the temp directory'
         print *, '               attempting to run in current dir'
      end if

c   o create "parallel" joda string
      write(xjoda,*) 'xjoda -procs ',nprocs,' -rank ',irank
#ifdef _DEBUG
      print *, '@PARENV_START: xjoda is "',xjoda,'"'
#endif

#else /* SERIAL BUILD */

      call aces_com_parallel_aces
      xjoda = 'xjoda'

#endif

      return
      end

