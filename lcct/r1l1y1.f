
      subroutine r1l1y1(icore,maxcor,uhf,listl1,listl1off,listr1,
     &                  listr1off)
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)
      common/sym/pop,vrt,nt,nfmi,nfea
C Input variables
      integer maxcor,uhf,listr1,listr1off,listl1,listl1off
C Pre-allocated local variables
      integer icore(maxcor)
C Local variables
      integer iofflist,imode,spin

      iofflist = 200

      imode = 0
      do spin = 1,uhf+1
C Create lists
        call inipck(1,8+spin,8+spin,53+spin+iofflist,imode,0,1)
        call inipck(1,8+spin,11-spin,55+spin+iofflist,imode,0,1)
        call inipck(1,10+spin,10+spin,57+spin+iofflist,imode,0,1)
        call updmoi(1,nfmi(spin),spin,91+iofflist,0,0)
        call updmoi(1,nfea(spin),spin,92+iofflist,0,0)
C Zero out the lists
        call zerolist(icore,maxcor,53+spin+iofflist)
        call zerolist(icore,maxcor,55+spin+iofflist)
        call zerolist(icore,maxcor,57+spin+iofflist)
        call zero(icore,max(nfmi(spin),nfea(spin)))
        call putlst(icore,1,1,1,spin,91+iofflist)
        call putlst(icore,1,1,1,spin,92+iofflist)
      end do

      call t1ring_r(icore,maxcor,uhf,.false.,iofflist,listr1,listr1off)
      call fiapart(icore,maxcor,uhf,listr1,listr1off,iofflist)
      call t1inw2_r(icore,maxcor,uhf,listr1,listr1off,iofflist)
      call f1inl1_r(icore,maxcor,uhf,iofflist,listl1,listl1off)
      call f2inl1_r(icore,maxcor,uhf,iofflist,listl1,listl1off)

      call l1inl1_r(icore,maxcor,uhf,1,iofflist,listl1,listl1off)
      if (uhf .ne. 0)
     &       call l1inl1_r(icore,maxcor,uhf,2,iofflist,listl1,listl1off)

      return
      end
