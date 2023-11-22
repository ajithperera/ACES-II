      subroutine stretch(density,nso,nao)

      implicit none
#include <aces.h>

      integer
     &   nso, nao, i, offset
      M_REAL
     &   density(nao*nao*2)

      offset=nao**2-nso**2
      do i=(nao**2+nso**2),(nao**2+1),(-1)
         density(i)=density(i-offset)
      end do

      return
      end
