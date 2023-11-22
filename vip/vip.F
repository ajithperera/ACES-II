#include "flags.h"
      Program vip
      Implicit none
      Integer iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"

      Call Aces_init(icore, i0,icrsiz, iuhf, .true.)
C
      Call Drive_vip(Icore(i0), icrsiz, Iuhf)
      
      Call Aces_fin
C
      Stop
      End


