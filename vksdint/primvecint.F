      subroutine primvecint

c This routine determines the occupancy and the MO to primitive function
c transformation matrix.  This is done in three steps:
c   1)  Get the orbitals (either SCF or natural) and determine occupancy
c   2)  Read the MO to non-symmetry-adapted AO transformation matrix
c       (written by the SCF code).
c   3)  Contract this with the primitive to non-symmetry-adapted AO
c       transformation matrix.

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <flags.com>
#include <kscore.com>
#include <mol.com>
#include <energ.com>
#include <intgrtflags.com>

#include <int.com>

      integer zcoeff,zscr,zscr2,zscr3,zscr4,i,
     &    zcmp2c1,pordr,setptr,znull,pnull

c zcoeff    : C1 MO coefficients
c pordr     : Scratch array for sorting MO vectors
c zcmp2c1   : transformation matrix from computational to c1 vectors
c zscr      : holds natural orbitals
c zscr2     : holds natural orbitals transformed to c1 AO's
c zscr3     : holds scf eigenvalues for reordering of correlated vectors
c zscr4     : holds occupation numbers

#ifdef _DEBUG_INTGRT
      write(*,*) '@INTGRT: Entered ', __FILE__
#endif /* _DEBUG_INTGRT */

      call callstack_push('PRIMVEC')

c Global data
      zxocc=setptr(1,F_REAL,2*nbastot)
      do i=1,2
        zpcoeff(i)=setptr(1,F_REAL,nbastot*totprim)
      end do

c Local data
      znull  =setptr(1,F_REAL,    1)
      pnull  =setptr(1,F_INTEGER, 1)

      zcoeff =setptr(1,F_REAL,    nbastot*naobasfn)
      zscr   =setptr(1,F_REAL,    nbastot*naobasfn)
      zscr2  =setptr(1,F_REAL,    nbastot*naobasfn)
      zscr3  =setptr(1,F_REAL,    nbastot)
      zscr4  =setptr(1,F_REAL,    nbastot)
      zcmp2c1=setptr(1,F_REAL,    nbastot*naobasfn)
      pordr  =setptr(1,F_INTEGER, naobasfn)

      if (memknown.eq.0) goto 999

      if (idns.eq.0) then
        if (.not.int_ks) write(*,9000)
      else
        if (iflags(18).eq.1) then
          ehar=echar
          if (.not.int_ks) write(*,9010)
        endif
      endif

 9000 format(/'Results using the SCF density')
 9010 format(/'Results using the correlated relaxed density')

c Skip this integration for correlated calculations with no relaxed
c density
      if (idns.eq.1.and.iflags(18).eq.0) goto 999

      do i=1,iuhf+1
        call gtaoAG(i,dcore(zscr),dcore(zscr2),
     &      dcore(zcmp2c1),kscore(pordr),dcore(zscr3),
     &      dcore(zcoeff),kscore(pnshellao),kscore(pnangatom),
     &      dcore(zscr4),dcore(zxocc))
        call xgemm('N','N',totprim,nbastot,naobasfn,1.d0,
     &      dcore(zprimcoef),totprim,dcore(zcoeff),naobasfn,0.d0,
     &      dcore(zpcoeff(i)),totprim)
      end do


c      call mkpra(dcore(zpcoeff(1)))

  999 continue
      call relptr(1,F_REAL,   znull)
cap      call relptr(1,F_INTEGER,pnull)

      call callstack_pop
      return
      end
