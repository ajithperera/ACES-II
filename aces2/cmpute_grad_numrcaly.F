
      subroutine cmpute_grad_numrcaly

#include "aces2.h"
#include "parallel_aces.com"

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

      call energy
      call c_gtod(is,ius)
           print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                        is+1.d-6*ius-dTimeStart,' seconds'
      return
      end

