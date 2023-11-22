#include "flags.h"

      program tdee
      implicit none

      integer iuhf
      Logical MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,CC2

      External Tdee_oed_pes

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

      Common /METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC

c ----------------------------------------------------------------------
      Call aces_init(icore, i0,icrsiz, iuhf, .true.)
C
      Call Tdee_driver(Icore(i0), icrsiz/iintfp, Iuhf)

      Call aces_fin
C
c ----------------------------------------------------------------------
      stop
      end

