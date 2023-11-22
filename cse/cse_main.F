#include "flags.h"

      program  cse
      implicit none

      integer iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
      Call aces_init(icore, i0,icrsiz, iuhf, .true.)

      Call cse_control(Icore(i0),icrsiz/iintfp,Iuhf)

      Call aces_fin
C
c ----------------------------------------------------------------------
      Stop
      End

