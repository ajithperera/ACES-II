
      subroutine FormSig(icore,maxcor,uhf,listr1,listr1off,listr2,
     &                   listr2rs)
      implicit none
C     Common Blocks
#include "machsp.com"
#include "flags.h"
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2)
      common/sym/pop,vrt,nt
      integer iflags(100)
      common/flags/iflags
      integer list2i1,list2i2,list2i3
      common /t2ilis/list2i1,list2i2,list2i3
C     Input Variables
      integer uhf,maxcor,listr1,listr1off,listr2,listr2rs
C     Pre-allocated Local Variables
      integer icore(maxcor)
C     External variables
      double precision dnrm2
      external dnrm2
C     Local Variables
      integer i000,i010,i020,kmax,mxcor,nbasis,iside,totlen,imode
      double precision signorm
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iside=1
      kmax=iflags(h_IFLAGS_estate_maxc)-1
C     Memory Pointers
      i000=1
      call getrec(20,'JOBARC','NBASTOT ',1,nbasis)

C     Load the relevant records and put them into a single list
C     for lineqy to utilize
      call maklst(icore(i000),maxcor,uhf)

C     These list numbers should be right.
      mxcor=maxcor/iintfp
      call loadvec2(1,icore(i000),mxcor,uhf,list2i3+1,0,list2i1-1,
     &              totlen,.false.)
      call dscal(totlen,-1.0d0,icore(i000),1)
      call aces_list_touch(1,totlen,1,373,0)
      call putlst(icore(i000),1,1,1,1,373)
      i010=i000+max(nbasis,kmax*kmax)*iintfp
      i020=i010+kmax*iintfp
      if (i020.gt.maxcor) call insmem('formsig',i020,maxcor)
      mxcor=(maxcor-i020+1)/iintfp
      call lineqy(icore(i020),mxcor,uhf,icore(i000),icore(i010),kmax,1,
     &            1,iside)

C Set up lists for Sigma
      imode = 0
      call inipck(1,13,14,listr2+2,imode,0,1)
      call inipck(1,9,9,listr2rs,imode,0,1)
      call inipck(1,10,9,listr2rs+2,imode,0,1)
      call inipck(1,9,10,listr2rs+3,imode,0,1)
      call inipck(1,11,12,listr2rs+5,imode,0,1)
      call updmoi(1,nt(1),1+listr1off,listr1,0,0)
      call aces_list_memset(1+listr1off,listr1,0)
      if (uhf .ne. 0) then
        call inipck(1,1,3,listr2,imode,0,1)
        call inipck(1,2,4,listr2+1,imode,0,1)
        call inipck(1,12,11,listr2rs+4,imode,0,1)
        call inipck(1,10,10,listr2rs+1,imode,0,1)
        call updmoi(1,nt(2),2+listr1off,listr1,0,0)
        call aces_list_memset(2+listr1off,listr1,0)
      endif
C Move Sigma amplitudes into lists listr2 and listr1
      call getlst(icore(i000),1,1,1,1,373)
      signorm = dnrm2(totlen,icore(i000),1)
      write(6,3000) signorm
      call updates(1,icore(i000),listr2,listr1off,listr1,uhf)
C Resort Sigma amplitudes
      call resort(icore(i000),maxcor,uhf,1,listr2,listr2rs)

      return
 3000 FORMAT(T3,'@FORMSIG, Norm of sigma amplitudes: ',F15.10,'.')
      end
