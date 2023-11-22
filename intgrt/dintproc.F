
      subroutine dintproc(ksdint,density,v,gradient,
     &                    nao,nso,
     &                    dscr,iscr,
     &                    order)

      implicit none
#include <aces.h>
#include <machsp.com>
#include <iuhf.com>
#include <mol.com>

      integer
     &   nao, nso,
     &   order(nao), iscr(nao)

      M_REAL
     &   ksdint(nao,nao,3,iuhf+1),
     &   density(nao,nao,iuhf+1),
     &   v(nao,nao),
     &   dscr(nao,nao),
     &   gradient(natoms*3)

      integer
     &   n, xyz, i, j, offset, spin

      M_REAL  ddot
      real sdot

      call callstack_push('DINTPROC')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c 1) Get the gradient.
c
c 2) Since we need the trace of the product PC'VC, transform the density
c    from SO/comp ordering to AO/zmat ordering once instead of
c    transforming every potential (of which there are 3*natoms).
c
c 3) Reorder ksdint from AO/angular to AO/zmat with the map KSSRTORD
c    which was written to JOBARC earlier.
c 
c 4) Separate the centers by looping over all atoms and extract from
c    ksdint only the elements that contribute to the particular
c    atom (loop index). Observe:
c
c       ksdint(mu,nu,q,s) = < mu | Vxc,s | d(nu)/dq >
c
c       v(mu,nu) =   < d(mu)/dq | Vxc,s | nu >
c                  + < mu | Vxc,s | d(nu)/dq >
c 
c 4) Add the dot product to the appropriate gradient element.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the gradient
      call getrec(1,'JOBARC','GRADIENT',natoms*3*iintfp,gradient)

#ifdef _DEBUG_VKSDINT
      write(*,*)
      write(*,*) "There are ",nao," AOs and ",nso," SOs."
      write(*,*)
      write(*,*) 'The old (computationally ordered) gradient:'
      write(*,*) 'atm    x               y               z'
 100  format     (i3 ,3f16.10)
      do n=1,natoms
         offset=3*n
         write(6,100) n, gradient(offset-2),
     &                   gradient(offset-1),
     &                   gradient(offset)
      end do
#endif /* _DEBUG_VKSDINT */

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the density matrix
      call mkdens(dscr,density)

c   o if this is a UHF calculation and if # SOs is less than # AOs then
c     we need to stretch out the matrix.
      if ((iuhf.eq.1).and.(nso.lt.nao)) call stretch(density,nso,nao)

c   o get the SO/comp to AO/zmat map
      call getrec(1,'JOBARC','CMP2ZMAT',nao*nso*iintfp,dscr)

c   o create CPC' using v (currently unused) as tmp
      do spin=1,iuhf+1
         call B_SYMM('R','L',nao,nso,
     &               1.d0,density(1,1,spin),nso,
     &                    dscr,             nao,
     &               0.d0,v,                nao)
         call B_GEMM('N','T',nao,nao,nso,
     &               1.d0,v,                nao,
     &                    dscr,             nao,
     &               0.d0,density(1,1,spin),nao)
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the AO/zmat to AO/ang map
      call getrec(1,'JOBARC','KSSRTORD',nao,order)

c   o reverse from AO/ang to AO/zmat ordering
      do spin=1,iuhf+1
         do xyz=1,3
            call mat_reorder
     &         (2,1,nao,nao,ksdint(1,1,xyz,spin),dscr,iscr,order)
         end do
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the AO center assignments in zmat ordering
      call getrec(1,'JOBARC','CNTERBF0',nao,order)

c   o get the computational to ZMAT ordering map
      call getrec(1,'JOBARC','MAP2ZMAT',natoms,iscr)

#ifdef _DEBUG_VKSDINT
      write(*,*) '@DINTPROC: iuhf is ',iuhf
      write(*,*) '@DINTPROC: CNTERBF0 is ',(order(n),n=1,nao)
      write(*,*) '@DINTPROC: MAP2ZMAT is ',(iscr(n),n=1,natoms)
      write(*,*) '@DINTPROC: nao is ',nao
      write(*,*) '@DINTPROC: density is ',density
      write(*,*) '@DINTPROC: ksdint is ',ksdint
#endif /* _DEBUG_VKSDINT */

      do n=1,natoms
         do spin=1,iuhf+1
            do xyz=1,3

c            o build v
               do j=1,nao
                  if (order(j).eq.iscr(n)) then
                     do i=1,nao
                        v(i,j)=ksdint(i,j,xyz,spin)
                     end do
                  else
                     do i=1,nao
                        v(i,j)=0.d0
                     end do
                  end if
               end do

c            o add the dot product to the gradient
               offset=3*n-3+xyz
               if (iuhf.eq.0) then
                  gradient(offset)=
     &            gradient(offset)
     &            +4*B_DOT(nao*nao,density(1,1,spin),1,v,1)
               else
                  gradient(offset)=
     &            gradient(offset)
     &            +2*B_DOT(nao*nao,density(1,1,spin),1,v,1)
               end if
#ifdef _DEBUG_VKSDINT
               write(*,*) '@DINTPROC: gradient(',offset,') += ',
     &            B_DOT(nao*nao,density(1,1,spin),1,v,1),'*2*(iuhf+1)'
#endif /* _DEBUG_VKSDINT */
            end do
         end do
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef _DEBUG_VKSDINT
      write(6,*) 'The new gradient:'
      write(6,*)
     & 'atm    x               y               z'
      do n=1,natoms
         offset=3*n
         write(6,100) n, gradient(offset-2),
     &                   gradient(offset-1),
     &                   gradient(offset)
      end do
#endif /* _DEBUG_VKSDINT */

c   o What are we waiting for?
      call putrec(1,'JOBARC','GRADIENT',natoms*3*iintfp,gradient)

      call callstack_pop
      return
      end

