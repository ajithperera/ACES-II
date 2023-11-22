#include "flags.h"

      Subroutine Drive_vip(Icore, Maxcor, Iuhf)
      
      Implicit Double Precision (A-H, O-Z)

      Dimension ICore(Maxcor)

#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "istart.com"

#ifdef _DEBUG_LVL0
      CALL AMPSUM(Icore,MAXCOR,IUHF,0,.TRUE.,'T')
#endif
C
      CALL INIT_IPEOM(Icore,Maxcor,Iuhf,Icalc)

      CALL SOLVE_IPEOM(Icore,Maxcor,Iuhf)

      RETURN
      END

