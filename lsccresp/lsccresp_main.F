#include "flags.h"

      program lsccresp
      implicit none

      integer iuhf,J0,Jcrsiz 
      integer Jcore(1)

      External Tdee_oed_pes

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
      Call aces_init(icore, i0,icrsiz, iuhf, .true.)

      Jcrsiz = 10000
      Jcore(1) = 0
      Call aces_malloc(icrsiz,Jcore,J0)
      if (Jcore(1).eq.0) then
         print *, '@-lsccresp_main: Unable to allocate at least ',
     &            Jcrsiz, ' integers of memory'
         call aces_exit(1)
      end if
C
      Call lsccresp_driver(Icore(i0),Jcore(J0),icrsiz/iintfp,Jcrsiz,
     &                     Iuhf)

      Call aces_fin
C
c ----------------------------------------------------------------------
      stop
      end

