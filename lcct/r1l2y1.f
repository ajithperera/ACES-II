
      subroutine r1l2y1(icore,maxcor,uhf,listl2,listr1,listr1off)
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nfmi(2),nfea(2)
      common/sym/pop,vrt,nt,nfmi,nfea
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer uhf,listr1,listr1off,maxcor,listl2
C Pre-allocated local variables
      integer icore(maxcor)
C Local variables
      integer iofflist,i0y(2),i0f(2),mxcor,i000,spin,imode
      double precision one
      data one /1.0d0/

C
C FOR CCSD SECOND DERIVATIVE CALCULATIONS, STORE W(IFMN) ON LISTS 207 TO 210
C
      iofflist = 200
      imode = 0
      call inipck(1,14,18,10+iofflist,imode,0,1)
      if (uhf .ne. 0) then
        call inipck(1,3,16,7+iofflist,imode,0,1)
        call inipck(1,4,17,8+iofflist,imode,0,1)
        call inipck(1,14,11,9+iofflist,imode,0,1)
      endif
      call inipck(1,13,11,30+iofflist,imode,0,1)
      if (uhf .ne. 0) then
        call inipck(1,1,9,27+iofflist,imode,0,1)
        call inipck(1,2,10,28+iofflist,imode,0,1)
        call inipck(1,13,18,29+iofflist,imode,0,1)
      endif

      call zerolist(icore,maxcor,10+iofflist)
      call zerolist(icore,maxcor,30+iofflist)
      if (uhf .ne. 0) then
        call zerolist(icore,maxcor,7+iofflist)
        call zerolist(icore,maxcor,8+iofflist)
        call zerolist(icore,maxcor,9+iofflist)
        call zerolist(icore,maxcor,27+iofflist)
        call zerolist(icore,maxcor,28+iofflist)
        call zerolist(icore,maxcor,29+iofflist)
      endif

c Form the W~ terms
      call sortring(icore,maxcor,uhf,1,1,123,124,118,117,125,126,.TRUE.)
      call dhbiajk3_r(icore,maxcor,uhf,1,1,listr1,listr1off,123,124,118,
     &                117,125,126,6+iofflist,.FALSE.)
      call sortring(icore,maxcor,uhf,1,2,123,124,118,117,125,126,.TRUE.)
C Forms Z(IA,JK) = R(AM) * W(IM,JK)
      call dhbiajk4_r(icore,maxcor,uhf,1,1,listr1,listr1off,50,
     &                6+iofflist)

C Expects FULL W(AB,CD)
      call w5ring_r(icore,maxcor,uhf,iofflist,listr1,listr1off)

      if (uhf .ne. 0) then
        call w5ab2_r(icore,maxcor,uhf,.false.,.true.,.false.,.false.,
     &               iofflist,listr1,listr1off)
        call w5ab1_r(icore,maxcor,uhf,.true.,.false.,.false.,.false.,
     &               iofflist,listr1,listr1off)
        call w5aa1_r(icore,maxcor,uhf,.true.,.false.,.false.,.false.,
     &               iofflist,listr1,listr1off)
      else
        call w5ab2rhf_r(icore,maxcor,uhf,.false.,.true.,.false.,.false.,
     &                  iofflist,listr1,listr1off)
      endif

C Now have formed the integrals to run L2 -> Y1

      i0y(1)=1
      if (uhf .ne. 0) then
        i0y(2) = i0y(1)+iintfp*nt(1)
      else
        i0y(2) = i0y(1)
      endif
      i000 = i0y(2)+iintfp*nt(2)
      mxcor = maxcor - i000 + 1
      do spin = 1,uhf+1
        call zero(icore(i0y(spin)),nt(spin))
      end do
C Contract L2 with the W~ intermediates

      if (uhf .eq. 1) then
        call t2t1aa2_r(icore(i0y(1)),icore(i000),mxcor,pop(1,1),
     &                 vrt(1,1),1,listl2,iofflist)
        call t2t1aa2_r(icore(i0y(2)),icore(i000),mxcor,pop(1,2),
     &                 vrt(1,2),2,listl2,iofflist)
      endif
      call t2t1ab2_r(icore(i0y(1)),icore(i000),mxcor,pop(1,1),pop(1,2),
     &               vrt(1,1),vrt(1,2),1,uhf,listl2,iofflist)
      
      if (uhf .eq. 1)
     &  call t2t1ab2_r(icore(i0y(2)),icore(i000),mxcor,pop(1,2),
     &                 pop(1,1),vrt(1,2),vrt(1,1),2,uhf,listl2,iofflist)

      if (uhf .eq. 1) then
        call t2t1aa1_r(icore(i0y(1)),icore(i000),mxcor,pop(1,1),
     &                 vrt(1,1),1,listl2,iofflist)
        call t2t1aa1_r(icore(i0y(2)),icore(i000),mxcor,pop(1,2),
     &                 vrt(1,2),2,listl2,iofflist)
      endif
      call t2t1ab1_r(icore(i0y(1)),icore(i000),mxcor,pop(1,1),pop(1,2),
     &             vrt(1,1),vrt(1,2),1,uhf,listl2,iofflist)
      if (uhf .eq. 1)
     &  call t2t1ab1_r(icore(i0y(2)),icore(i000),mxcor,pop(1,2),
     &                 pop(1,1),vrt(1,2),vrt(1,1),2,uhf,listl2,iofflist)

C Three-body terms
      call gformg(1,1,listl2,44,100,icore(i000),mxcor,0,1.0d0,uhf)
      i0f(1) = i000
      i0f(2) = i0f(1)+uhf*iintfp*nt(1)
      i000 = i0f(2) + iintfp*nt(2)
      call getlst(icore(i0f(1)),1,1,1,1,93+iofflist)
      if (uhf .ne. 0) call getlst(icore(i0f(2)),1,1,1,2,93+iofflist)
      mxcor = maxcor - i000 + 1
      call r1l2y1a(icore(i0y(1)),icore(i0f(1)),icore(i000),mxcor,uhf)
      call r1l2y1b(icore(i0y(1)),icore(i0f(1)),icore(i000),mxcor,uhf)

      do spin = 1,uhf+1
        call getlst(icore(i0f(1)),1,1,1,2+spin,90)
        call daxpy(nt(spin),1.0d0,icore(i0f(1)),1,icore(i0y(spin)),1)
        call putlst(icore(i0y(spin)),1,1,1,2+spin,90)
      end do
      return
      end
