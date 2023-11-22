      subroutine setupgrd(max_angpts)

c This subroutine sets up the radial grid and the angular grid for each
c angular point.

      implicit none
#include <aces.h>
#include <kscore.com>
#include <intgrtflags.com>
#include <mol.com>
#include <grid.com>

      integer igrid,setptr

      integer max_angpts

c The routine(s) initGRIDTYPE(gridnum,subtype) (where GRIDTYPE=lebed or
c some other grid type) should do the following:
c    if the first argument is 0
c       set maxangpt
c           maxanggrd
c           ntotrad
c    else
c       set grdangpts(igrid)
c           gridxyz(3,maxangpts,igrid)
c           gridwt(maxangpts,igrid)
c    end if

#ifdef _DEBUG_INTGRT
      write(6,*) '@INTGRT: Entered ', __FILE__
#endif /* _DEBUG_INTGRT */
          memknown =-1

      call callstack_push('SETUPGRD')

c When determining the amount of memory needed, we have to initialize
c the radial and angular grids before we set the pointers to determine
c their sizes.
      if (memknown.ne.1) then
        numangfct=0
        call initgrid(kscore)
          write(6,*) 'what sis up'
        numgrid=1
        maxangpts=0
        numgrid=0
        call rdradgrd(kscore)
        ntotrad=0
        do igrid=1,numgrid
           write(6,*) gridlist(igrid,1)
          if (gridlist(igrid,1).eq.1) then
            call initlebed(0,igrid,gridlist(igrid,2),kscore,dcore,dcore)
          endif
        end do
      endif
      pradgrid  =setptr(1,F_INTEGER, int_numradpts)
      pgrdangpts=setptr(1,F_INTEGER, numgrid)
      zgridxyz  =setptr(1,F_REAL,    3*maxangpts*numgrid)
      zgridwt   =setptr(1,F_REAL,    maxangpts*numgrid)
      pangfct   =setptr(1,F_INTEGER, 3*numangfct)
     
c The second time through the program, we need to set the pointers before
c we read in the full description of the radial and angular grids.
        write(6,*) 'jjjj'
      if (memknown.ne.0) then
          write(6,*) ' what is inside'
        call initgrid(kscore(pangfct))
        call rdradgrd(kscore(pradgrid))
        do igrid=1,numgrid
          if (gridlist(igrid,1).eq.1) then
            call initlebed(1,igrid,gridlist(igrid,2),kscore(pgrdangpts),
     &          dcore(zgridxyz),dcore(zgridwt))
          endif
        end do
      endif

c Set up an ordered list of symmetry unique centers (polist) and
c which axes in the integration grid the axes in the Aces2 coordinate
c system correspond to (ix, iy and iz)

      polist =setptr(1,F_INTEGER, compnorb)
      pintegaxis = setptr(1,F_INTEGER, natoms*3)

      if (memknown.ne.0) then
        call ordlst(kscore(polist),kscore(pcompmemb),kscore(pcomppopv),
     &      dcore(zcoord),kscore(pintegaxis))
      endif

c Determine factors that will stay constant throughout the entire
c integration.

      zatmvc =setptr(1,F_REAL, natoms*natoms*3)
      zrij   =setptr(1,F_REAL, natoms*natoms)
      zatmpth=setptr(1,F_REAL, natoms*natoms*3*minpt)
      zptdis =setptr(1,F_REAL, natoms*natoms*minpt)
      zprsqrd=setptr(1,F_REAL, natoms*natoms*minpt*natoms)
      zpthpt =setptr(1,F_REAL, natoms*natoms*3*minpt*natoms)
      zbslrd =setptr(1,F_REAL, natoms)
      zaij   =setptr(1,F_REAL, natoms*natoms)

      if(memknown.ne.0) then
        call intconst(kscore(patomchrg),dcore(zcoord),dcore(zatmvc),
     &      dcore(zrij),dcore(zatmpth),dcore(zptdis),dcore(zprsqrd),
     &      dcore(zpthpt),dcore(zbslrd),dcore(zaij))
      endif

      max_angpts=maxangpts
      call callstack_pop
      return
      end
