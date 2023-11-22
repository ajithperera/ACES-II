
      program main
      implicit none

#include "icore.com"
#include "istart.com"
#include "method.com"

      call crapsi(icore,iuhf,0)
      call nmr(icore(i0),icrsiz,iuhf)
      call aces_fin

      stop
      end

