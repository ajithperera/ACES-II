#include "flags.h"
      program a2mix_main
      implicit none
      Integer iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"

      Call aces_init(icore, i0,icrsiz, iuhf, .true.)
C
      Call a2mix(Icore(i0), icrsiz, Iuhf)

      Call aces_fin
C
      Stop
      End


