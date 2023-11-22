      subroutine fixv(v,ord,screx,scr)

c On entry v conatins the matrix elements from EXX.
c Add the elements coming from the numerical integration and
c find the SA AO Fock matrix containing exchange-correlation
c potential(s).

#include <aces.h>
      implicit none
#include <iuhf.com>
#include <machsp.com>
#include <mol.com>

      integer ord(naobasfn,2),ispin

      M_REAL
     &    v(naobasfn,naobasfn,iuhf+1),
     &    screx(naobasfn,naobasfn),
     &    scr(naobasfn,naobasfn,iuhf+1)

      integer nnsa,nsa,i,j

      call callstack_push('FIXV')

      nnsa=naobasfn
      nsa=nbastot

      call getrec(1,'JOBARC','KSSRTORD',nnsa,ord(1,1))
c      call getrec(1,'JOBARC','KSPOTENT',nnsa*nnsa*(iuhf+1)*iintfp,scr)

      do ispin=1,iuhf+1
           do i=1,nnsa-1
           do j=i+1,nnsa
             scr(i,j,ispin)=scr(j,i,ispin)
           end do
         end do
c      call dprt(nnsa,nnsa,scr(1,1,ispin),'V(nsa,nsa)')

       call mat_reorder(2,1,nnsa,nnsa,scr(1,1,ispin),screx(1,1),
     &   ord(1,2),ord(1,1))

c      call dcopy(nnsa*nnsa,scr(1,1,ispin),1,v(1,1,ispin),1) 

c       do i=1,nnsa
c         do j=1,nnsa
c            v(i,j,ispin)=v(i,j,ispin)+scr(i,j,ispin)
c         end do
c       end do

      end do

      call daxpy(nnsa*nnsa*(iuhf+1),1.d0,scr,1,v,1)
      call getrec(1,'JOBARC','CMP2ZMAT',nnsa*nsa*iintfp,screx(1,1))

c V(sa) = Yt x V(nsa) x Y

      do ispin=1,iuhf+1
        call B_GEMM('t','n',nsa,nnsa,nnsa,1.d0,screx(1,1),nnsa,
     &    v(1,1,ispin),nnsa,0.d0,scr(1,1,2),nnsa)
        call B_GEMM('n','n',nsa,nsa,nnsa,1.d0,scr(1,1,2),nnsa,
     &    screx(1,1),nnsa,0.d0,v(1,1,ispin),nsa)

c        call dprt(nsa,nsa,v(1,1,ispin),'V(sa,sa)')
      end do

      call callstack_pop
      return
      end
