
      program main
      implicit none

#include "icore.com"
#include "istart.com"
#include "machsp.com"
      integer iuhf

      call crapsi(icore,iuhf,0)
      call props(icore(i0),icrsiz/iintfp,icore(i0),icrsiz)
      call aces_fin

      stop
      end

