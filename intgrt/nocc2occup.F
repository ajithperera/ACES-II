            subroutine nocc2occup

C   this routine will occupy nocc to occup for dft works

#include<aces.h>
         implicit none
#include<iuhf.com>
#include<symm2_ks.com>
         integer i,jj
         integer NOCC(16)

         common /POPUL/NOCC



         do i=1,8
             occup(i,1)=NOCC(i)
         end do
 
          do jj=1,8
             occup(jj,2)=NOCC(8+jj)
          end do
         
         return
         end




