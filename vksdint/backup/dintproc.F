      subroutine dintproc(ksdint,v,density,gradient,
     &                    nao,nso,
     &                    dscr,iscr,
     &                    order)

#include <aces.h>
#include <io.par>
#include <callstack.com>
#include <const.par>
#include <machsp.com>
#include <intgrtflags.com>
#include <flags.com>
#include <mem.com>
#include <mol.com>
#include <energ.com>
#include <grid.com>
#include <int.com>
#include <dftfunc.com>

      M_REAL
     &   ksdint(nao,nao,3,iuhf+1),
     &   density(nao,nao,iuhf+1),
     &   v(nao,nao),
     &   dscr(nao,nao),
     &   gradient(natoms*3)

      integer
     &   nao, nso,
     &   order(nao), iscr(nao)

      integer
     &   n, xyz, i, j, offset, spin

      M_REAL  ddot
      integer sdot

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
      call getrec(1,'GRADIENT',natoms*3*iintfp,gradient)

#ifdef DEBUGACES
      write(6,*) ''
      write(6,*) 'The old (computationally ordered) gradient:'
      write(6,*) 'atm    x               y               z'
 100  format     (i3 ,3f16.10)
      do n=1,natoms
         offset=3*n
         write(6,100) n, gradient(offset-2),
     &                   gradient(offset-1),
     &                   gradient(offset)
      end do
#endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the density matrix
      call mkdens(dscr,density)

      write(*,*) "There are ",nao," AOs and ",nso," SOs."
      do spin=1,iuhf+1
      do j=1,nso
      do i=1,nso
         write(*,*) "Pso(",i,",",j,",",spin,"): ",density(i,j,spin)
      end do
      end do
      end do

c   o get the SO/comp to AO/zmat map
      call getrec(1,'CMP2ZMAT',nao*nso*iintfp,dscr)

c   o create CPC' using v (currently unused) as tmp
      do spin=1,iuhf+1
         call B_SYMM('R','L',nao,nso,
     &               one, density(1,1,spin),nso,
     &                    dscr,             nao,
     &               zero,v,                nao)
         call B_GEMM('N','T',nao,nao,nso,
     &               one, v,                nao,
     &                    dscr,             nao,
     &               zero,density(1,1,spin),nao)
         do j=1,nao
         do i=1,nao
            write(*,*) "Pao(",i,",",j,",",spin,"): ",density(i,j,spin)
         end do
         end do
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the AO/zmat to AO/ang map
      call getrec(1,'KSSRTORD',nao,order)

c   o reverse from AO/ang to AO/zmat ordering
      do spin=1,iuhf+1
         do xyz=1,3
            call mat_reorder
     &         (2,1,nao,nao,ksdint(1,1,xyz,spin),dscr,iscr,order)
         end do
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c   o get the AO center assignments in zmat ordering
      call getrec(1,'CNTERBF0',nao,order)

c   o get the computational to ZMAT ordering map
      call getrec(1,'MAP2ZMAT',natoms,iscr)

      do n=1,natoms
         do spin=1,iuhf+1
            do xyz=1,3

C DOLT! Since we sum over all mu and nu and since P is symmetric, we
C only need to take 2*P*ksdint instead of P*(ksdint(i,j)+ksdint(j,i)).
C BEGIN OLD PIECE
Cc            o build v
C               do j=1,nao
C
Cc               o the < (mu) | Vxc | d(nu)/d(XYZ) > piece
C                  if (order(j).eq.iscr(n)) then
C                     do i=1,nao
C                        v(i,j)=ksdint(j,i,xyz,spin)
C                     end do
C                  else
C                     do i=1,nao
C                        v(i,j)=zero
C                     end do
C                  end if
C
Cc               o the < d(mu)/d(XYZ) | Vxc | (nu) > piece
C                  do i=1,nao
C                     if (order(i).eq.iscr(n))
C     &                  v(i,j)=v(i,j)+ksdint(i,j,xyz,spin)
C                  end do
C               end do
C END OLD PIECE
C BEGIN NEW PIECE
c            o build v
               do j=1,nao
                  if (order(j).eq.iscr(n)) then
                     do i=1,nao
                        v(i,j)=ksdint(i,j,xyz,spin)
                     end do
                  else
                     do i=1,nao
                        v(i,j)=zero
                     end do
                  end if
               end do
C END NEW PIECE

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
            end do
         end do
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef DEBUGACES
      write(6,*) 'The new gradient:'
      write(6,*)
     & 'atm    x               y               z'
      do n=1,natoms
         offset=3*n
         write(6,100) n, gradient(offset-2),
     &                   gradient(offset-1),
     &                   gradient(offset)
      end do
#endif

c   o What are we waiting for?
      call putrec('GRADIENT',natoms*3*iintfp,gradient)

      call callstack_pop
      return
      end

