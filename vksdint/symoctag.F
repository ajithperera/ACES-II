      subroutine symoctag
     &    (atomchrg,atmvc,rij,aij,cdnt,rsqrd,rrtmp,wtintr,
     &    totwt,radgrid,radpt,rwt,integaxis,gridxyz,gridwt,grdangpts,
     &    iradpt,iangpt,grid,evalpt,dw,dp,dz,ppp,dmmuji)

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
     &    rrtmp(natoms),radpt(int_numradpts),
     &    rwt(int_numradpts,maxangpts),
     &    gridxyz(3,maxangpts,numgrid),gridwt(maxangpts,numgrid)
      logical evalpt

      integer 
     &    i,n,ip,io,in

      M_REAL
     &    atmwt,fctr,xx(0:3),totwt

      M_REAL
     & dw(3,natoms),dp(3,natoms,natoms),dz(3,natoms),ppp(natoms),
     & dmmuji(3)


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
c      if (xx(integaxis(1,icntr)).ge.0.d0.and.
c     &    xx(integaxis(2,icntr)).ge.0.d0.and.
c     &    xx(integaxis(3,icntr)).ge.0.d0) then
        
        call crtgridAG(icntr,atomchrg,atmvc,rij,aij,cdnt,rsqrd,
     &      rrtmp,wtintr,atmwt,xx(1),gridwt,rwt,iangpt,iradpt,grid,
     &      dw,dp,dz,ppp,dmmuji)

c Determine total weight
        totwt=
c     &       pi*4.d0*atmwt*rwt(iradpt,iangpt)*
c     &      gridwt(iangpt,grid)
     &        atmwt*rwt(iradpt,iangpt)*gridwt(iangpt,grid)


c       if ( iangpt .eq. 1 ) then
c       write(*,*)iradpt,totwt,atmwt,rwt(iradpt,iangpt),gridwt(iangpt,1)
c       write(*,*) 'icnt=',icntr
c       end if
c Multiply weight by number of symmetry equivalent points
c        n=0
c        if (xx(integaxis(1,icntr)).ne.0.d0) n=n+1
c        if (xx(integaxis(2,icntr)).ne.0.d0) n=n+1
c        if (xx(integaxis(3,icntr)).ne.0.d0) n=n+1
c        fctr=2.d0**n

         fctr=1.d0
c        write(*,*) 'fctr=',fctr,totwt
        totwt=fctr*
     &  totwt 

            
c      else
c        evalpt=.false.
c      endif
     
      call callstack_pop
      return
      end
