      subroutine symoct
     &    (atomchrg,atmvc,rij,aij,cdnt,rsqrd,rrtmp,wtintr,
     &    totwt,radgrid,radpt,rwt,integaxis,gridxyz,gridwt,grdangpts,
     &    iradpt,iangpt,grid,evalpt)

c This routine sets up the cartesian coordinates for the numerical 
c integration grid

      implicit none
#include <aces.h>
#include <intgrtflags.com>
      M_REAL pi
#include <mol.com>
#include <grid.com>
#include <int.com>
      integer
     &    atomchrg(natoms),radgrid(int_numradpts),
     &    integaxis(3,natoms),grdangpts(numgrid),
     &    iradpt,iangpt,grid,ncount
      M_REAL
     &    cdnt(natoms,3),atmvc(natoms,natoms,3),
     &    rsqrd(natoms),rij(natoms,natoms),
     &    aij(natoms,natoms),wtintr(natoms),
c     &    itotwt(icntr,int_numradpts,maxangpts),
     &    rrtmp(natoms),radpt(int_numradpts),
     &    rwt(int_numradpts,maxangpts),
     &    gridxyz(3,maxangpts,numgrid),gridwt(maxangpts,numgrid)
      logical evalpt

      integer 
     &    i,n,ip,io,in

      M_REAL
     &    atmwt,fctr,xx(0:3),totwt

c ********************************************************************
c ********************************************************************

      call callstack_push('SYMOCT')
      pi = acos(-1.d0)

c Determine the angular weights and the cartesian coordiantes of each point
      xx(0)=0.d0

c Determine the cartesian coordinates of the point with respect to the
c integration center
      do 185 i=1,3
        xx(i)=radpt(iradpt)*gridxyz(i,iangpt,grid)
  185 continue

c Only keep points which fall in the xyz octant
      if (xx(integaxis(1,icntr)).ge.0.d0.and.
     &    xx(integaxis(2,icntr)).ge.0.d0.and.
     &    xx(integaxis(3,icntr)).ge.0.d0) then
        
        call crtgrid(icntr,atomchrg,atmvc,rij,aij,cdnt,rsqrd,
     &      rrtmp,wtintr,atmwt,xx(1))

c Determine total weight
        totwt=
     &       pi*4.d0*atmwt*rwt(iradpt,iangpt)*
     &      gridwt(iangpt,grid)

c Multiply weight by number of symmetry equivalent points
        n=0
        if (xx(integaxis(1,icntr)).ne.0.d0) n=n+1
        if (xx(integaxis(2,icntr)).ne.0.d0) n=n+1
        if (xx(integaxis(3,icntr)).ne.0.d0) n=n+1
        fctr=2.d0**n

        totwt=fctr*
     &  totwt 

      else
        evalpt=.false.
      endif
      
      call callstack_pop
      return
      end
