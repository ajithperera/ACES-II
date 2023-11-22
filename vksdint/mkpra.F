       subroutine mkpra(CA)

       implicit none
#include <aces.h>
#include <iuhf.com>
#include <mol.com>

        M_REAL
     &  CA(totprim,naobasfn)

         integer i,j

CSSS       write(*,*) totprim
CSSS       write(*,*) naobasfn
CSSS       do i=1,totprim
CSSS         write(*,*)(CA(i,j),j=1,naobasfn),i
CSSS       end do
      

      return
      end
 
