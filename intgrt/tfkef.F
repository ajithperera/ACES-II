      function tfkef(roa,rob,grdaa,grdab,grdbb)
      
c Determine the value of the Thomas Fermi kinetic energy functional.
      
      implicit none
#include <aces.h>
#include <callstack.com>

      M_REAL
     &    roa,rob,grdaa,grdab,grdbb

      M_REAL
     &    tfkef,ro,tfcnst,pi

      callstack_curr='TFKEF'
      pi = acos(-1.d0)
      ro=roa+rob
      tfcnst=(3.d0/10.d+00)*(3.d0*pi*pi)**(2.d0/3.d0)
      tfkef=tfcnst*(ro**(5.d0/3.d0))
     
      return
      end
