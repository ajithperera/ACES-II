#include "flags.h"
                  
      program main
      implicit none

c INTERNAL VARIABLES
      integer iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"

c ----------------------------------------------------------------------

      call start_time
      call aces_init_rte
      call aces_ja_init
      call getrec(1,'JOBARC','IFLAGS', 100,iflags)
      call getrec(1,'JOBARC','IFLAGS2',500,iflags2)

c ----------------------------------------------------------------------

      icrsiz = iflags(h_IFLAGS_mem)
      icore(1) = 0
      do while ((icore(1).eq.0).and.(icrsiz.gt.1000000))
         call aces_malloc(icrsiz,icore,i0)
         if (icore(1).eq.0) icrsiz = icrsiz - 1000000
      end do

c----------------------------------------------------------------------
      call aces_io_init(icore,i0,icrsiz,.true.)
      if (iflags(h_IFLAGS_ref) .eq. 0) then    
         iuhf = 0
      else       
        iuhf = 1
      endif     
c
      call aces_init_chemsys
C
      call TrIIII(icore(i0), icrsiz, iuhf)
C
      call aces_fin
      Stop
      End
