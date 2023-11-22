      function fuzzy(n,f0)

#include <aces.h>
#include <const.par>
#include <callstack.com>

      integer n,i
      M_REAL
     &    fuzzy,f0,f1

      callstack_curr='FUZZY'
      do i=1,n
        f1=three*half*f0-half*f0**three
        f0=f1
      end do
      fuzzy=f0

      return
      end
