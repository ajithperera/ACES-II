      function fuzzy(n,f0)

      implicit none
#include <aces.h>
#include <callstack.com>

      integer n,i
      M_REAL
     &    fuzzy,f0,f1

      callstack_curr='FUZZY'
      do i=1,n
        f1=3.d0*0.5d0*f0-0.5d0*f0*f0*f0
        f0=f1
      end do
      fuzzy=f0

      return
      end
