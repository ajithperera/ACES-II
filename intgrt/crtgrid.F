      subroutine crtgrid(icntr,atomchrg,atmvc,rij,aij,cdnt,rsqrd,
     &    rrtmp,wtintr,atmwt,xx)

c This routine determines the cartesian coordinates for arbitrary grid
c point

      implicit none
#include <aces.h>
#include <intgrtflags.com>
#include <mol.com>
#include <grid.com>
#include <callstack.com>
      integer icntr,atomchrg(natoms)
      M_REAL
     &    cdnt(natoms,3),atmvc(natoms,natoms,3),
     &    rsqrd(natoms),rij(natoms,natoms),
     &    aij(natoms,natoms),wtintr(natoms),
     &    rrtmp(natoms),atmwt,xx(3)

      integer i,iatom,jatom
      M_REAL
     &    wttot,zmuij,xmuij,fuzzy,cutij,f4

      callstack_curr='CRTGRID'

c Determine the cartesian coordinates of the point with respect to
c every atom
      do 180 iatom=1,natoms
        rsqrd(iatom)=0.d0
        do 170 i=1,3
          cdnt(iatom,i)=atmvc(icntr,iatom,i)+xx(i)
          rsqrd(iatom)=rsqrd(iatom)+cdnt(iatom,i)**2
  170   continue

        rrtmp(iatom)=dsqrt(rsqrd(iatom))

  180 continue

      if(int_parttyp.eq.rigid)then
        atmwt=1.d0
      else
        wttot=0.d0

        do 200 iatom=1,natoms
          if (atomchrg(iatom).eq.110) then
            wtintr(iatom)=0.d0
            goto 200
          else
            wtintr(iatom)=1.d0
          endif

          do 190 jatom=1,natoms
            if (atomchrg(jatom).eq.110) goto 190
            if (iatom.eq.jatom) goto 190

            zmuij=(rrtmp(iatom)-rrtmp(jatom))/
     &          rij(iatom,jatom)
c Hideo found this bug from Becke's paper
c                 xmuij=zmuij +aij(iatom,jatom)*(1.d0-aij(iatom,jatom)*
c     &                aij(iatom,jatom))
            xmuij=zmuij +aij(iatom,jatom)*(1.d0-zmuij**2)
            f4=fuzzy(int_fuzzyiter,xmuij)
            cutij=0.5d0*(1.d0-f4)
            wtintr(iatom)=wtintr(iatom)*cutij

c trial ***********************************************************
c                wtintr(iatom)=(0.5d0-0.5d0*zmuij)*wtintr(iatom)
c *****************************************************************

  190     continue
          wttot=wttot+wtintr(iatom)
  200   continue
        atmwt=wtintr(icntr)/wttot

      endif

      return
      end
