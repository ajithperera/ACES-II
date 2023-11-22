      subroutine aux_trans(nlda,transf,ineigv,orig,res,scrau,
     &                     order,scrorder)

#include <aces.h>
#include <const.par>
#include <machsp.com> 
#include <kscore.com>
#include <mol.com>

      integer nlda,order(naobasfn),scrorder(naobasfn)

      M_REAL
     &    transf(nlda,nlda),
     &    ineigv(nlda),
     &    orig((naobasfn*(naobasfn+1))/2,nlda),
     &    res((naobasfn*(naobasfn+1))/2,nlda),
     &    scrau(naobasfn,naobasfn)

      integer i,j,k,l,naonao

      call callstack_push('AUX_TRANS')

      naonao=(naobasfn*(naobasfn+1))/2

      open (unit=46,file='OVERINX',status='OLD',access='SEQUENTIAL',
     &    form='UNFORMATTED')
       rewind(46)
       read (46) orig
       close(46,status='DELETE')
c

      call xgemm('n','n',naonao,nlda,nlda,
     &    one,orig,naonao,transf,nlda,zero,
     &    res,naonao)

c
c     .. debugging the KSPLOT module                                  ..
c
c$$$      call mat_print (
c$$$     $      orig, naonao, nlda, '3 index overlap after the orthog.',
c$$$     $      'ZAO', 'LDA', 6, 3)
c
c     .. To get rid of linearly dependent aux. basis functions, the   ..
c     .. auxilliary basis set has been orthogonolized.        Let's   ..
c     .. save the transformation matrix in case we want to plot the   ..
c     .. OEP local potential later                                    ..
c
c     .. attention here! we also need to zero out some columns if     ..
c     .. there are linearly depedent functions. NOT IMPLEMENTED YET   ..
c
              
      open (
     $     unit=44, file='AUXTRANSF',
     $     status='new', form='unformatted')
      write (44) ((transf (i,j), i=1,nlda), j=1, nlda)
      close (44, status='keep')
c
      open (unit=47,file='KSAOINT',status='new',
     &    form='unformatted')
          rewind(47)

      do k=1,nlda
        if (ineigv(k).ne.zero) then
          i=1
          j=1
          do l=1,naonao
            scrau(i,j)=res(l,k)
            scrau(j,i)=scrau(i,j)
            j=j+1
            if (j.gt.naobasfn) then
               i=i+1
               j=i
            end if
          end do

          call mat_reorder(2,1,naobasfn,naobasfn,scrau,orig,
     &         scrorder,order)
          do i=1,naobasfn
            do j=1,naobasfn
               write(47) scrau(i,j)
            end do
          end do
        end if
      end do
      close(47,status='KEEP')

      call callstack_pop
      return
      end
