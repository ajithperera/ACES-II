      program nlorb
      implicit none

      integer iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"

      call aces_init(icore,i0,icrsiz,iuhf,.true.)
      call drive(icore(i0),icrsiz,iuhf)

      call aces_fin

      stop
      end
