
      subroutine final_check

      implicit none
#include "aces.h"
#include "mol.com"

#ifdef _DEBUG_INTGRT
      write(*,*) '@FINAL_CHECK: totprim is ',totprim
#endif /* _DEBUG_INTGRT */

      return
      end

