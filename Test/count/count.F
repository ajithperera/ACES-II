
      program main
      implicit none

c INTERNAL VARIABLES
c   o the effective dcore size
      integer dcoredim

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"

c ----------------------------------------------------------------------

c   o initialize xcc and allocate memory
      call xcc_init(dcoredim)

c   o pass the address to a driver so we don't have to worry about
c     scaling values by iintfp
      call xcc_driver(dcoredim,icore(i0))

      call crapso
      stop 0
      end

