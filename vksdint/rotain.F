         subroutine rotain(gradient,coord)

         implicit none
#include <aces.h>
#include <machsp.com>
#include <iuhf.com>
#include <mol.com>

          M_REAL
     &  gradient(natoms,3),coord(3,natoms)

      M_REAL
     &  rotx,roty,rotz
      integer i
      M_REAL
     &  transx,transy,transz
 
       rotx=0.d0
       roty=0.d0
       rotz=0.d0
       transx=0.d0
       transy=0.d0
       transz=0.d0

c       write(*,*) 'coord=',coord(1,1),coord(1,2)

       do i=1,natoms
       rotx=rotx+coord(1,i)*gradient(i,2)
     &        -coord(2,i)*gradient(i,1)
       roty=roty+coord(2,i)*gradient(i,3)-
     &         coord(3,i)*gradient(i,2)
       rotz=rotz+coord(3,i)*gradient(i,1)
     &        -coord(1,i)*gradient(i,3)
       transx=transx+gradient(i,1)
       transy= transy+gradient(i,2)
       transz=transz+gradient(i,3)
       end do

CSSS       write(*,*) 'rotin=',rotx,roty,rotz
CSSS       write(*,*)'tran=',transx,transx,transz
      end

