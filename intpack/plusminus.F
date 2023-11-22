      subroutine plusminus(iexp,wt,x,y,z,grid,gridxyz,gridwt)

c This routine takes a point (x,y,z) with it's weight (wt) and adds it (in
c all possible negative and positive combinations) to the expanded list of
c grid points.  The point(s) should be added starting at the iexp'th location
c in the list.

#include <aces.h>
#include <const.par>
#include <grid.com>
#include <callstack.com>

      integer iexp,grid
      M_REAL
     &    wt,x,y,z,
     &    gridxyz(3,maxangpts,numgrid),gridwt(maxangpts,numgrid)

      integer ix,iy,iz,mx,my,mz

      callstack_curr='PLUSMINUS'

c We'll perform do loops from 1 to 1 or -1 to 1 (step 2) to get all
c possible combinations.
      mx=-1
      if (x.eq.zero) mx=1
      my=-1
      if (y.eq.zero) my=1
      mz=-1
      if (z.eq.zero) mz=1

      do 30 ix=mx,1,2
        do 20 iy=my,1,2
          do 10 iz=mz,1,2

            gridwt(iexp,grid)=wt
            gridxyz(1,iexp,grid)=ix*x
            gridxyz(2,iexp,grid)=iy*y
            gridxyz(3,iexp,grid)=iz*z
            iexp=iexp+1

   10     continue
   20   continue
   30 continue

      return
      end
