      subroutine initgrid(angfct)

c This initializes some constants used in the grids.

#include <aces.h>
#include <grid.com>
#include <maxangshell.par>
      integer angfct(numangfct,3)

      integer i,j,k,x,y,z,ptr

c Atomic size using Slater's rules for the radial integration
      data atmrad /1.000,0.588,
     &    3.077,2.051,1.538,1.231,1.026,0.879,0.769,0.684,
     &    4.091,3.158,2.571,2.169,1.875,1.651,1.475,1.333,
     &    6.27,4.84,4.59,4.38,4.20,4.01,3.82,3.68,3.53,
     &    1.192,3.27,3.16,2.76,2.44,2.19,1.98,1.81,1.66,
     &    7.29,5.67,5.37,5.12,4.91,4.69,4.46,4.31,4.12,
     &    3.97,3.82,3.68,3.23,2.85,2.55,2.31,2.12,1.94,
     &    8.03,6.24,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     &    0.00,0.00,0.00,0.00,0.00,0.00,0.00,5.90,5.60,
     &    5.40,5.16,4.91,4.74,4.54,4.37,4.20,4.06,3.55,
     &    3.14,2.80,0.00,0.00,0.00/

c Bragg-Slater radii for determining the relative size of the
c polyhedra in the polyatomic integration scheme
      data xbsl /0.35,0.50,
     &    1.45,1.05,0.85,0.70,0.65,0.60,0.50,0.65,
     &    1.80,1.50,1.25,1.10,1.00,1.00,1.00,0.95,
     &    2.20,1.80,1.60,1.40,1.35,1.40,1.40,1.40,1.35,
     &    0.60,0.60,1.35,1.30,1.25,1.15,1.15,1.15,1.10,
     &    2.35,2.00,1.80,1.55,1.45,1.45,1.35,1.30,1.35,
     &    1.40,1.60,1.55,1.55,1.45,1.45,1.40,1.40,1.30,
     &    2.60,2.15,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
     &    0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.75,1.55,
     &    1.45,1.35,1.35,1.30,1.35,1.35,1.35,1.50,1.90,
     &    1.75,1.60,0.00,0.00,0.00/

      call callstack_push('INITGRID')

c Length of angular momentum function vectors
      write(6,*) 'numangfct=',numangfct,maxangshell
      if (numangfct.eq.0) then
        do i=1,maxangshell
          numangfct=numangfct+i*(i+1)/2
        end do
        goto 999
      end if

c The number of angular momentum functions
c angfct(*,1) = 0  100  211000  3221110000  ...
c angfct(*,2) = 0  010  010210  0102103210  ...
c angfct(*,3) = 0  001  001012  0010120123  ...
      ptr=1
      do i=1,maxangshell
        x=i-1
        do j=1,i
          y=j-1
          z=0
          do k=1,j
            angfct(ptr,1)=x
            angfct(ptr,2)=y
            angfct(ptr,3)=z
            y=y-1
            z=z+1
            ptr=ptr+1
          end do
          x=x-1
        end do
      end do

  999 continue
      call callstack_pop
      return
      end
