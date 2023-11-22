C add_ecp_int.F: communicates between CSD and ECP code
C                normal integrals computed

      subroutine add_ecp_int(VNA,coords)

      use ecpcom

      implicit none

      include "lengths.h"
      include "CI1PRB.h"
      include "CI2PRB.h"

C In/out parameters
      double precision :: VNA(norbp,nat), coords(3,nat)
      integer :: i, j

C      print*,"VNA before"
C      print*,VNA
C      print*
C      print*,"Max VNA before",maxval(VNA)
C      print*,"Min VNA before",minval(VNA)

C     get the normal integrals
      call ecp_glue(coords,VNA,norb,nat)

C      print*,"VNA after"
C      print*,VNA
C      print*
C      print*,"Max VNA after",maxval(VNA)
C      print*,"Min VNA after",minval(VNA)

      return
      end subroutine

