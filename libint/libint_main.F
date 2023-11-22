#include "flags.h"

      Program Libint_main
      implicit None

      integer Iuhf,J0,Jcrsiz
      integer Jcore(1)

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
      Call Aces_init(Icore,I0,Icrsiz,Iuhf,.True.)

      Jcrsiz = 10000000
      Jcore(1) = 0
      Call aces_malloc(icrsiz,Jcore,J0)
      if (Jcore(1).eq.0) then
         print *, '@Libint_main: Unable to allocate at least ',
     &            Jcrsiz, ' integers of memory'
         call aces_exit(1)
      end if
C
      Call Libint_driver(Icore(i0),Jcore(J0),icrsiz/Iintfp,
     &                   Jcrsiz/Iintfp,Iuhf)

      Call Aces_fin
C
c ----------------------------------------------------------------------
      stop
      end

