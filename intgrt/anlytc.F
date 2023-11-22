      subroutine anlytc

c This routine drives determining analytic values for the kinetic
c energy, nuclear attraction energy, coulomb energy and exchange
c energy using the HF and correlated density.

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <vmol.com>
#include <flags.com>
#include <intgrtflags.com>
#include <kscore.com>
#include <mol.com>
#include <int.com>
#include <energ.com>
      integer setptr

      integer zdens,zcoef,zbuf,pbuf,zdens1,zcdens,pmap,zocc,
     &    znull,pnull,length 
      M_REAL xx

#ifdef _DEBUG_INTGRT
      write(*,*) '@INTGRT: Entered ', __FILE__
#endif /* _DEBUG_INTGRT */

      int_printlev=1
      if (int_ks) then
        if (int_ks_finaliter) then
          int_printlev=2
        else
          int_printlev=3
        endif
      endif

      call callstack_push('ANLYTC')

c zdens       Alpha/Beta SCF density
c zcoef       Scratch for holding SCF MO vectors
c zbuf        Scratch for holding AO integrals
c zdens1      Alpha/beta density matrix contracted for one electron integrals
c zcdens      Alpha/beta density matrix contracted for coulomb integrals
c zocc        Natural orbital occupation

c pbuf        Scratch for holding AO integral labels
c pmap        Array for mapping C1 to symmetry

      znull  =setptr(1,F_REAL,    1)
      pnull  =setptr(1,F_INTEGER, 1)

      zcoef  =setptr(1,F_REAL,    nbastot*nbastot)
      zdens  =setptr(1,F_REAL,    nbastot*nbastot*2)
      zbuf   =setptr(1,F_REAL,    ilnbuf)
      zdens1 =setptr(1,F_REAL,    nbastot*(nbastot+1)/2*(iuhf+1))
      zcdens =setptr(1,F_REAL,    nbastot*nbastot*(iuhf+1))
      zocc   =setptr(1,F_REAL,    nbastot)

      pbuf   =setptr(1,F_INTEGER, ilnbuf)
      pmap   =setptr(1,F_INTEGER, nbastot)

      if (memknown.eq.0) goto 999
      call mkdens(dcore(zcoef),dcore(zdens))
      if (idns .eq. 0) call exchange(dcore(zbuf),kscore(pbuf),
     &    dcore(zdens),
     &    dcore(zdens1),dcore(zcdens),
     &    kscore(pmap),ehfkin,ehfnatr,ehfx,
     &    ehfcoul,ehar,escf,0)
C#ifdef _NOSKIP
C I don't think this block is need
      call getrec(0,'JOBARC',"AONTORBA",length,kscore)
      if (iflags(2).ne.0 .and. iflags(2).ne.41 .and.
     &    iflags(18).eq.1 .and. length .gt.0) then
        call mknatden(dcore(zcoef),dcore(zocc),dcore(zdens))
        call exchange(dcore(zbuf),kscore(pbuf),dcore(zdens),
     &      dcore(zdens1),
     &      dcore(zcdens),kscore(pmap),
     &      eckin,ecnatr,ecx,eccoul,echar,ecor,1)

CSSS        write(*,9100)
CSSS        xx=ecnatr-ehfnatr
CSSS        write(*,9110) xx
CSSS        xx=eckin-ehfkin
CSSS        write(*,9120) xx
CSSS        xx=eccoul-ehfcoul
CSSS        write(*,9130) xx
CSSS        xx=ecx-ehfx
CSSS        write(*,9140) xx
CSSS        xx=ecor-escf
CSSS        write(*,9150) xx

      endif
C#endif 

 9100 format(/'The correlation contribution ')
 9110 format('   Nuclear-electron attraction energy   =',f20.12)
 9120 format('   Kinetic energy                       =',f20.12)
 9130 format('   Coulomb energy                       =',f20.12)
 9140 format('   Exchange energy                      =',f20.12)
 9150 format('   Total Hartree-Fock energy expression =',f20.12/)

  999 continue
      call relptr(1,F_REAL,   znull)
cap      call relptr(1,F_INTEGER,pnull)

      call callstack_pop
      return
      end
