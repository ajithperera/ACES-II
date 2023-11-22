
      subroutine form_dens_lamsig(icore,maxcor,uhf,nonhf,listl1,
     &                            listl1off,listl2,listl2rs,listr1,
     &                            listr1off,listr2,listr2rs,listgrl,
     &                            listgtl,listgrlof,listgtlof,listtmp,
     &                            listtmpof,listt1,listt1off)
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer maxcor,uhf,listl1,listl1off,listl2,listl2rs,listr1,
     &        listr1off,listr2,listr2rs,listgrl,listgtl,listgrlof,
     &        listgtlof,listtmp,listtmpof,listt1,listt1off
      logical nonhf
C Preallocated local variables
      integer icore(maxcor)
C Local variables
      integer idoo,idvv,idvo,idov,itop,mxcor,ioffoo,ioffvv,ioffvo,
     &        ioffov,spin,imode
      double precision znorm,fact

      znorm = 0.0d0
      idoo = 1
      idvv = idoo + (nd1(1)+uhf*nd1(2))*iintfp
      idvo = idvv + (nd2(1)+uhf*nd2(2))*iintfp
      idov = idvo + (nt(1)+uhf*nt(2))*iintfp
      itop = idov + (nt(1)+uhf*nt(2))*iintfp
      mxcor = maxcor - itop

      call updmoi(1,nt(1),1+listtmpof,listtmp,0,0)
      call aces_list_memset(1+listtmpof,listtmp,0)
      if (uhf .ne. 0) then
        call updmoi(1,nt(2),2+listtmpof,listtmp,0,0)
        call aces_list_memset(2+listtmpof,listtmp,0)
      endif
C Form one-particle density matrix and update the densities
      call gdens(1,1,1,icore(idoo),icore(idvv),icore(idvo),icore(idov),
     &           icore(itop),mxcor,uhf,znorm,listgrl,listgtl,listgrlof,
     &           listgtlof,listtmp,listtmpof,listr1,listl1,listr1off,
     &           listl1off,listr2,listl2,listr2rs,listl2rs,listt1,
     &           listt1off,nonhf)
C Update densities on disk
      ioffoo = idoo
      ioffvv = idvv
      ioffvo = idvo
      ioffov = idov
      fact = 1.0d0
      do spin = 1,uhf+1
        call updens(icore(ioffoo),icore(ioffvv),icore(ioffvo),
     &              icore(ioffov),icore(itop),spin,fact,nonhf)
        ioffoo = ioffoo+nd1(1)*iintfp
        ioffvv = ioffvv+nd2(1)*iintfp
        ioffvo = ioffvo+nt(1)*iintfp
        ioffov = ioffov+nt(1)*iintfp
      end do
C Form two-particle density matrix for <0|(1+Lambda)Hbar Sigma|0>
      call inigam(uhf,icore,maxcor)

      call tpdens_lcc(icore,maxcor,uhf,listr1,listl1,listr1off,
     &                listl1off,listr2,listl2,listr2rs,listl2rs)

      return
      end
