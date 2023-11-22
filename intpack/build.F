      subroutine build(iordr,nshellao,nangmom)

#include <aces.h>
#include <io.par>
#include <const.par>
#include <flags.com>
#include <intgrtflags.com>
#include <machsp.com>
#include <sym.com>
#include <mol.com>
#include <int.com>
#include <kscore.com>

      integer
     &    iordr(naobasfn),nangmom(natoms),
     &    nshellao(maxshell,natoms)

      integer
     &    iangmom,ifrom,
     &    istart,iatom,shell,nao,nshell,ito,iao,ishell

      call callstack_push('BUILD')

c Put higher angular momentum vectors in correct order
c Coming in, the p orbitals are ordered (x,y,z,x,y,z,...)
c After reordering, they are ordered (x,x,...,y,y,...,z,z,...)

      istart=1
      do iatom=1,natoms
        shell=1
        do iangmom=1,nangmom(iatom)
          nao=nshellao(shell,iatom)
          if (nao.eq.0) goto 30
          nshell=iangmom*(iangmom+1)/2
          shell=shell+nshell
          ito=istart
          do iao=1,nao
            ifrom=istart+iao-1
            do ishell=1,nshell
              iordr(ito)=ifrom
              ifrom=ifrom+nao
              ito=ito+1
            end do
          end do
          istart=istart+nao*nshell
   30     continue
        end do
      end do

      call callstack_pop
      return
      end
