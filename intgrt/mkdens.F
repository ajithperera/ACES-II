      subroutine mkdens(coef,dens)

c This routine creates the SCF density needed to determine the
c SCF kinetic energy, nuclear attraction energy, total one 
c electron energy, coulomb energy, hartree energy and exchange
c energy.

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <machsp.com>
#include <flags.com>
#include <symm2_ks.com>
#include <mol.com>
      M_REAL
     &    coef(nbastot,nbastot),dens(nbastot,nbastot,iuhf+1)

      integer ioff1,ioff2,irrep,nocct,iocc,i,j,ispin
      character*1 type(2)
      data type/'A','B'/

      call callstack_push('MKDENS')

      do 50 ispin=1,iuhf+1

c Get the matrix of mo's
        call getrec(1,'JOBARC','SCFEVC'//spinc(ispin)//'0',
     &      iintfp*nbastot*nbastot,coef)

c Total density
        ioff1=0
        ioff2=0
        do 40 irrep=1,nirrep
          nocct=occup(irrep,ispin)
          do 30 iocc=1,nocct
            do 20 i=1,numbasir(irrep)
              do 10 j=1,numbasir(irrep)
                dens(ioff2+j,ioff2+i,ispin)=
     &              dens(ioff2+j,ioff2+i,ispin)+
     &              coef(ioff2+i,ioff1+iocc)*coef(ioff2+j,ioff1+iocc)
   10         continue
   20       continue
   30     continue

c Problems...  SI
c 
c          if(iflags(2).eq.0.or.iflags(2).eq.41) then
            ioff1=ioff1+numbasir(irrep)
c          else
c            ioff1=ioff1+nocct
c          endif
          ioff2=ioff2+numbasir(irrep)
   40   continue

   50 continue

      call callstack_pop
      return
      end
