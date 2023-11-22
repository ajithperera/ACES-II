#include "f_types.h"
      Subroutine Getmem_a2(Icore,iCoreNdx,iCoreDim,iMem,bAllocMem)

      implicit none

#include "flags.com"
#include "machsp.com"
      integer              xFlag
      common /aces_reflag/ xFlag
      save   /aces_reflag/
c ARGUMENTS
      integer iCore(*), iCoreNdx, iCoreDim, iUHF
      logical bAllocMem

c EXTERNAL FUNCTIONS
      F_ADR c_adr
      external c_adr

c PARAMETERS
      integer iMemMin, iMemInc
      parameter (iMemMin=1*1024*1024,iMemInc=1*1024*1024)

c INTERNAL VARIABLES
      double precision dTmp
      integer i0, iMem, iTmp
 
c   o allocate memory and initialize the I/O subsystem
      if (bAllocMem) then

c      o allocate core memory
         if (xFlag.eq.0) then
            iCore(1) = 0
            do while ((iCore(1).eq.0).and.(iMem.gt.iMemMin))
               call aces_malloc(iMem,iCore,i0)
               if (iCore(1).eq.0) iMem = iMem - iMemInc
            end do
#ifdef _DEBUG_ACES_INIT
            print *, '@ACES_INIT: heap size   = ',iMem,' integers'
            print *, '            iCore index = ',i0
#endif
            if (iMem.lt.iflags(36)) then
               print *, '@ACES_INIT: MEMORY WARNING!'
               print *, '            requested ',iflags(36),' integers'
               print *, '            allocated ',iMem,' integers'
            end if
            if (iCore(1).ne.0) then
#ifdef _DEBUG
c            o make sure iCore(i0) starts on a double word boundary
               print *, '@ACES_INIT: heap range',c_adr(iCore(i0)),
     &                                           c_adr(iCore(i0+iMem))-1
#ifdef _PTRS_ARE_WORDS
               iTmp = iand(c_adr(iCore(i0)),iintfp-1)
               if (iTmp.ne.0) then
                  print *, '@ACES_INIT: alignment defect in aces_malloc'
                  print *, '            absolute offset = ',iTmp,' W'
                  call aces_exit(1)
               end if
#else
               iTmp = iand(c_adr(iCore(i0)),ifltln-1)
               if (iTmp.ne.0) then
                  print *, '@ACES_INIT: alignment defect in aces_malloc'
                  print *, '            absolute offset = ',iTmp,' B'
                  call aces_exit(1)
               end if
#endif
#endif
               xFlag = 1
            else
c               print *, '@ACES_INIT: Request for ',iMem,
c     &                  ' integers of memory failed.'
               print *, '@ACES_INIT: unable to allocate at least ',
     &                  iMemMin,' integers of memory'
               call aces_exit(1)
            end if
         end if

c      o initialize the I/O subsystem ('T' creates and initializes the
cache)
CSSS         call aces_io_init(iCore,i0,iMem,.true.)

c      o transfer the iCore statistics
         iCoreNdx = i0
         iCoreDim = iMem

c     else if (.not.bAllocMem)
      else

         iCore(1) = 0
         iCoreNdx = 1
         iCoreDim = 1

c     end if (bAllocMem)
      end if

      Return
      End
