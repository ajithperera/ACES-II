
c This routine reads the command-line arguments and initializes the
c parallel environment.

      subroutine xcc_cl_proc
      implicit none

#include "mpif.h"

c EXTERNAL FUNCTIONS
      integer f_strpbrk, c_atol
      integer*4 iargc

c INTERNAL VARIABLES
      character arg*80
      integer iLast
      integer*4 i, nArgs
      integer iMPIErr

c COMMON BLOCKS
#include "flags.com" /* for iflags(44) */
#include "xcc.com" /* for CNTOT and CNID */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ', __FILE__
#endif

c ----------------------------------------------------------------------

c   o assume a single process
      CNTOT = 1
      CNID  = 0
      call MPI_COMM_SIZE(MPI_COMM_WORLD,CNTOT,iMPIErr)
      if (iMPIErr.ne.0) then
         print *, '@XCC_CL_PROC: unable to retrieve MPI comm size'
         call xcc_errex
      end if
      call MPI_COMM_RANK(MPI_COMM_WORLD,CNID,iMPIErr)
      if (iMPIErr.ne.0) then
         print *, '@XCC_CL_PROC: unable to retrieve MPI comm rank'
         call xcc_errex
      end if

c   o override the number of cache slots
      iflags(44) = 100
      iflags(44) = 384

c ----------------------------------------------------------------------

c   o get the number of CL arguments
      nArgs = iargc()
      if (nArgs.lt.1) return

c   o loop over arguments
      do i = 1, nArgs
         call getarg(i,arg)

c      o CNTOT
         if (arg(1:2).eq.'-n') then
            call getarg(i+1,arg)
            iLast = f_strpbrk(arg,' ')
            arg(iLast:iLast) = achar(0)
            CNTOT = c_atol(arg)
         end if

c      o CNID
         if (arg(1:2).eq.'-i') then
            call getarg(i+1,arg)
            iLast = f_strpbrk(arg,' ')
            arg(iLast:iLast) = achar(0)
            CNID = c_atol(arg)
         end if

c      o cache slots
         if (arg(1:3).eq.'-cs') then
            call getarg(i+1,arg)
            iLast = f_strpbrk(arg,' ')
            arg(iLast:iLast) = achar(0)
            iflags(44) = c_atol(arg)
         end if

c   o end loop over arguments
c     end do i = 1, nArgs
      end do

c ----------------------------------------------------------------------

c VERIFY CONSISTENCY

c   o CNTOT
      if (CNTOT.lt.1) then
         print *, '@XCC_CL_PROC: resetting number of compute nodes to 1'
         CNTOT = 1
      end if

c   o CNID
      if (CNID.lt.0) then
         print *, '@XCC_CL_PROC: resetting compute node ID to 0'
         CNID = 0
      end if
      if (CNID.ge.CNTOT) then
         print *, '@XCC_CL_PROC: resetting compute node ID to ',CNTOT-1
         CNID = CNTOT-1
      end if

c   o cache slots
      if (iflags(44).lt.1) then
         print *, '@XCC_CL_PROC: resetting cache slots to 1'
         iflags(44) = 1
      end if

      print *, '@XCC_CL_PROC: summary'
      print *, '              (ID,TOTAL) = ',CNID,CNTOT
      print *, '              CACHE SLOTS = ',iflags(44)

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ', __FILE__
#endif

      return
      end

