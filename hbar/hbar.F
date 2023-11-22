#include "flags.h"

      program hbar
      implicit none

      integeriuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
      Call aces_init(icore, i0,icrsiz, iuhf, .true.)
C
      Call Drive_Hbar(Icore(i0), icrsiz/iintfp, Iuhf)

      Call aces_fin
C
c ----------------------------------------------------------------------
      stop
      end

