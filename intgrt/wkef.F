      function wkef(roa,rob,grdaa,grdab,grdbb)
      
c Determine the value of the Weizsacker kinetic energy functional.
      
      implicit none
#include <aces.h>
#include <callstack.com>

      M_REAL
     &    roa,rob,grdaa,grdab,grdbb

      M_REAL
     &    wkef,ro

      callstack_curr='WKEF'
      ro=roa+rob
      wkef=0.d0
      if (ro.ge.1.d-12) wkef=0.125d0*grdaa/ro

      return
      end
