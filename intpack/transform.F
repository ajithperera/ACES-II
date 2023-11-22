       subroutine transform(dimaubs,valint,tran,scr,scr2,eigv,
     &                      iordr,scrordr)

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

      integer dimaubs,iordr(naobasfn),scrordr(naobasfn)

      M_REAL
     &    valint((naobasfn*(naobasfn+1))/2,dimaubs),
     &    scr((naobasfn*(naobasfn+1))/2,dimaubs),
     &    scr2(naobasfn,naobasfn),
     &    tran(dimaubs,dimaubs),
     &    eigv(dimaubs)


      integer i,j,ldep,naonao,naux,nwork

      call callstack_push('TRANSFORM')

c Build the upper triangle for four-index intergrals
c or the overlap needed for three-index integrals.

      if (int_overlp.eq.4) then
         do i=1,dimaubs-1
            do j=i+1,dimaubs
               valint(i,j)=valint(j,i)
            end do
          end do
c          call dprt(dimaubs,dimaubs,valint,'INTEG')
      else
          do i=1,naobasfn-1
             do j=i+1,naobasfn  
                scr2(i,j)=scr2(j,i)
             end do
          end do
      end if

      open (unit=2,file='OVERINX',status='NEW',access='SEQUENTIAL',
     &    form='UNFORMATTED')
      rewind(2)
      write (2) valint
      close(2,status='KEEP')

      naonao=(naobasfn*(naobasfn+1))/2
      nwork=max(naonao*dimaubs, 4*dimaubs)
     

      if (int_overlp.eq.4) then
        call dgeev('N','V',dimaubs,valint,dimaubs,eigv,
     &     scr2,scr,naonao,tran,dimaubs,scr,nwork,j)

        call invert_vector(dimaubs,eigv,two,scr,ldep)
        do i=1,dimaubs
           eigv(i)=scr(i,i)
        end do

        call xgemm('n','n',dimaubs,dimaubs,dimaubs,one,
     &      tran,dimaubs,scr,dimaubs,zero,valint,naonao)

        call aux_trans(dimaubs,valint,eigv,tran,
     &                 scr,scr2,iordr,scrordr)
      else
        call dgeev('N','V',dimaubs,scr2,dimaubs,eigv,
     &     valint,scr,naonao,tran,dimaubs,scr,nwork,j)

        call invert_vector(dimaubs,eigv,two,scr2,ldep)
        do i=1,dimaubs
           eigv(i)=scr2(i,i)
           scr2(i,i)=zero
        end do

        do i=1,dimaubs
           do j=1,dimaubs
             scr2(j,i)=scr2(j,i)+tran(j,i)*eigv(i)
           end do
        end do

        call aux_trans(dimaubs,scr2,eigv,valint,
     &                 scr,tran,iordr,scrordr)   
      end if

      write(6,333) ldep

 333  format(/,'   There are ',i4,' linear dependencies in the Auxiliary
     & Basis Set') 
 
      naux=dimaubs-ldep
c      call putrec('NAUXBAS',1,       naux)   
      open(unit=3,file='NUMAUX',status='NEW',form='FORMATTED')
      write(3,999) naux
      close(3,status='KEEP')
 999  format(i5)
      call callstack_pop

      return
      end
