C add_ecp_deriv.F: communicates between CSD and ECP code
C                  calculates integral derivatives

      subroutine add_ecp_deriv(VNAD,coords)

      use ecpcom

      implicit none

      include "lengths.h"
      include "CI1PRB.h"
      include "CI2PRB.h"

C In/out parameters
      double precision :: VNAD(norbp,nat,6), coords(3,nat)

C      print*,"Max VNAD before",maxval(VNAD)
C      print*,"Min VNAD before",minval(VNAD)
      call ecp_glue_deriv(coords,VNAD,norb,nat)
C      print*,"Max VNAD before",maxval(VNAD)
C      print*,"Min VNAD before",minval(VNAD)

      return
      end subroutine

