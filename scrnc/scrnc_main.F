#include "flags.h"

      program scrnc 
      implicit none

      integer iuhf

      External Tdee_oed_pes

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
      Call aces_init(icore, i0,icrsiz, iuhf, .true.)
C
      Call scrnc_driver(Icore(i0), icrsiz/iintfp, Iuhf)

      Call aces_fin
C
c ----------------------------------------------------------------------
      stop
      end

