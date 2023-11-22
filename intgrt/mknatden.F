      subroutine mknatden(coef,occ,dens)

c This routine creates the correlated density needed to determine the
c correlated kinetic energy, nuclear attraction energy, total one
c electron energy, coulomb energy, hartree energy and exchange
c energy.

      implicit none
#include <aces.h>
#include <machsp.com>
#include <iuhf.com>
#include <mol.com>
      M_REAL
     &    coef(nbastot,nbastot),dens(nbastot,nbastot,iuhf+1),
     &    occ(nbastot)

      integer iocc,i,j,ispin

      call callstack_push('MKNATDEN')

      call dzero(dens,nbastot*nbastot*(iuhf+1))

      do 50 ispin=1,iuhf+1

c Get the matrix of mo's
        call getrec(1,'JOBARC','AONTORB'//spinc(ispin),
     &      iintfp*nbastot*nbastot,coef)
        call getrec(1,'JOBARC','OCCNUM_'//spinc(ispin),
     &      iintfp*nbastot,occ)

        if(iuhf.eq.0) then
          do 10 i=1,nbastot
            occ(i)=0.5d0*occ(i)
   10     continue
        endif

c Total density
        do 40 iocc=1,nbastot
          do 30 i=1,nbastot
            do 20 j=1,nbastot
              dens(j,i,ispin)=dens(j,i,ispin)+
     &            occ(iocc)*coef(i,iocc)*coef(j,iocc)
   20       continue
   30     continue
   40   continue

   50 continue

      call callstack_pop
      return
      end
