      subroutine numint

#include <aces.h>
#include <io.par>
#include <const.par>
#include <machsp.com>
#include <intgrtflags.com>
#include <flags.com>
#include <kscore.com>
#include <mol.com>
#include <grid.com>
#include <int.com>

      integer ncount,ncnt,nnn,ifctr,zirrtmp,zradpt,zrwt,
     &    zwtintr,zsize,zrelsiz,zrtmp,zrint,setptr,iradpt,iangpt,
     &    grid,zvalmo,zcdnt,zrsqrd,zvalprim,zxnat,
     &    znull,pnull,zvalao,
     &    zvalint,pordr,zscrord,zscr,zscr2,zeigval

       integer ztran

      M_REAL
     &    factor,totwt

      logical evalpt

      integer i,dimauxbs,naonao

      M_REAL
     &    times(20)

      call callstack_push('NUMINT')

      znull  =setptr(1,F_REAL,    1)
      pnull  =setptr(1,F_INTEGER, 1) 

      naonao=(naobasfn*(naobasfn+1))/2
      if (int_overlp.eq.4) then
         dimauxbs=(naobasfn*(naobasfn+1))/2
      else
         dimauxbs=naobasfn
      end if    

       write(6,*) 'naobasfn',naobasfn, dimauxbs
c zcdnt     : the coordinates of the current point with respect to each atom
c zrsqrd    : the radial distance sqared from each atom to the current point
c zvalprim  : the value of each primitive function at the current point
c zvalao    : the value of each AO at the current point
c zvalint   : the value of each ijkl overlap combination
c zwtintr   : An intermediate in calculating the weight due to each atom.
c zradpt    : The array of radial points

      zcdnt    =setptr(1,F_REAL, natoms*3)
      zrsqrd   =setptr(1,F_REAL, natoms)
      zvalprim =setptr(1,F_REAL, totprim)
      zxnat  =setptr(1,F_REAL, nbastot)
      zirrtmp=setptr(1,F_REAL, natoms)
      zradpt =setptr(1,F_REAL, int_numradpts)
      zrwt   =setptr(1,F_REAL, int_numradpts*maxangpts)
      zwtintr=setptr(1,F_REAL, natoms)
      zsize  =setptr(1,F_REAL, natoms*natoms)
      zrelsiz=setptr(1,F_REAL, natoms*natoms)
      zrtmp  =setptr(1,F_REAL, natoms)
      zrint  =setptr(1,F_REAL, natoms*maxangpts)

      zvalao =setptr(1,F_REAL, naobasfn)
      zvalint=setptr(1,F_REAL, naonao*dimauxbs)

      pordr  =setptr(1,F_INTEGER, naobasfn)
      zscrord=setptr(1,F_INTEGER, naobasfn)

      zscr   =setptr(1,F_REAL, naonao*dimauxbs)
      zscr2=setptr(1,F_REAL,naobasfn*naobasfn)
      zeigval=setptr(1, F_REAL,dimauxbs)
      ztran=setptr(1,F_REAL,dimauxbs*dimauxbs)

      if (memknown.eq.0) go to 999

      call dzero(dcore(zvalint),naonao*dimauxbs)

      if (int_overlp.eq.3) 
     &   call dzero(dcore(zscr2),naobasfn*naobasfn)

c Determine relative atomic size for the atomic partitioning
      call intpath(kscore(pnshellatom),kscore(pnshellprim),
     &    dcore(zalpha),dcore(zprsqrd),dcore(zvalprim),
     &    dcore(zpcoeff(1)),dcore(zvalmo),
     &    dcore(zpthpt),dcore(zsize),dcore(zptdis),
     &    dcore(zrelsiz),dcore(zbslrd),kscore(pangfct))

c Loop over integration centers (nuclei) and set the radial grid using
c Slater's rules.

c Set up ordered list of unique atoms
      if (c1symmet.eq.1) then
        ncount=natoms
      else
        ncount=fullnorb
      endif

      do i=1,20
        times(i)=zero
      end do

      do 30 ncnt=1,ncount
        nnn=ncnt

c Determine center and symmetry factor (number of symmetry equivalent atoms)
        if (c1symmet.eq.1) then
          icntr=ncnt
          ifctr=1
          factor=one
        else
          icntr=kscore(polist+ncnt-1)
          ifctr=kscore(pfullpopv+ncnt-1)
          factor=dble(ifctr)
        endif
c Set up grid for integrating only symmetry unique octants
        call oct(kscore(patomchrg),dcore(zatmvc),dcore(zrij),
     &      kscore(pradgrid),dcore(zrelsiz),dcore(zrtmp),dcore(zrint),
     &      dcore(zradpt),dcore(zrwt),dcore(zbslrd),
     &      dcore(zgridxyz),kscore(pgrdangpts))

        do iradpt=1,int_numradpts
          grid=kscore(pradgrid+iradpt-1)
          do iangpt=1,kscore(pgrdangpts+grid-1)
            evalpt=.true.

            call symoct(kscore(patomchrg),dcore(zatmvc),dcore(zrij),
     &          dcore(zaij),dcore(zcdnt),dcore(zrsqrd),dcore(zirrtmp),
     &          dcore(zwtintr),totwt,kscore(pradgrid),
     &          dcore(zradpt),dcore(zrwt),kscore(pintegaxis),
     &          dcore(zgridxyz),dcore(zgridwt),kscore(pgrdangpts),
     &          iradpt,iangpt,grid,evalpt)

c Integrate over center, icntr
            if (evalpt) then
                call integ(kscore(pnshellatom),kscore(pnshellprim),
     &          dcore(zcdnt),dcore(zalpha),dcore(zrsqrd),
     &          dcore(zvalprim),
     &          dcore(zpcoeff(1)),
     &          totwt,
     &          dcore(zxnat),dcore(zvalint),dcore(zscr),dcore(zscr2),
     &          factor,ifctr,
     &          dcore(zvalao),
     &          dcore(zprimcoef),dcore(zcoord),
     &          kscore(pfullmemb),
     &          kscore(pangfct),iradpt,iangpt,dimauxbs,
     &          times)
            end if
          end do
        end do

   30 continue

      call build(kscore(pordr),kscore(pnshellao),kscore(pnangatom))

      call transform(dimauxbs,dcore(zvalint),dcore(ztran),
     &     dcore(zscr),dcore(zscr2),dcore(zeigval),
     &     kscore(pordr),kscore(zscrord))

  999 continue

      call relptr(1,F_REAL,   znull)
      call relptr(1,F_INTEGER,pnull)

      call callstack_pop

      return
      end
