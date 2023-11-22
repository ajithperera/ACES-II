
      subroutine conv_check_wrap(icore,maxcor,uhf,conv,listv,
     &                           cycle)
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer maxcor,uhf,listv,cycle
C Output variables
      logical conv
C Pre-allocated local variables
      integer icore(maxcor)
C Local variables
      integer lenvo,ivai,iiter,maxiter,iend,ifai

      lenvo=max(nt(1),uhf*nt(2))
      ivai = 1
      iiter = ivai + lenvo*iintfp
      call getrec(20,'JOBARC','OOMBPTCY',1,maxiter)
      iend = iiter + maxiter*iintfp
      ifai = iend
      iend = ifai + lenvo*iintfp

      if (iend .gt. maxcor) call insmem('CONV_CHECK',iend,maxcor)

      call conv_check(icore(ivai),lenvo,icore(iiter),maxiter,cycle,uhf,
     &                conv,listv,icore(ifai))

      return
      end

      subroutine conv_check(vai,len,iterations,maxiter,cycle,uhf,conv,
     &           listv,fai)
      implicit none
#include "flags.h"
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer iflags(100)
      common/flags/iflags
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer uhf,len, maxiter,cycle,listv
      double precision vai(len),fai(len)
C Output variables
      logical conv
C Pre-allocated Local variables
      double precision iterations(maxiter)
C Local variables
      integer spin,ii
      double precision zilch,errmax,tol,fact,errmax1
      data zilch /0.0d0/

      tol = 10.0d0**(-dble(iflags(h_IFLAGS_cc_conv)))
C Check on the Frobenious norm

      errmax = zilch
      fact = 2.0d0 * dble(2-uhf)
      do spin = 1,uhf+1
        call getlst(vai,1,1,1,spin,listv)
        do ii = 1,nt(spin)
          errmax = errmax + fact*abs(vai(ii))**2
        end do
      end do
      errmax  = sqrt(errmax)

C If errmax < tol we are converged
      conv = errmax .le. tol
      call getrec(-1,'JOBARC','OOMBITER',(maxiter+1)*iintfp,iterations)
      iterations(cycle+1) = errmax
      call putrec(20,'JOBARC','OOMBITER',(maxiter+1)*iintfp,iterations)
      write(6,*) 'OO-MBPT Convergence Check: ',errmax

C Test Fai (only relevant on the first iteration, but whatever)
      errmax = zilch
      do spin = 1,uhf+1
        call getlst(vai,1,1,1,2+spin,93)
        do ii = 1,nt(spin)
          errmax = errmax + fact*abs(vai(ii))**2
        end do
      end do
      errmax = sqrt(errmax)
      write(6,*) 'Fai Frobenious Norm: ',errmax
      call getrec(-1,'JOBARC','OOMBITR2',(maxiter+1)*iintfp,iterations)
      iterations(cycle+1) = errmax
      call putrec(20,'JOBARC','OOMBITR2',(maxiter+1)*iintfp,iterations)

      return
      end
