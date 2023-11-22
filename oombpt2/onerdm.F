
      subroutine onerdm(icore,maxcor,uhf,singles)
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer maxcor,uhf
      logical singles
C Preallocated local variables
      integer icore(maxcor)
C Local variables
      integer id,iboff,iend,mxcor,itmp,spin

C OO Part of density
C All off-diagonal elements are zeroed out
      id = 1
      iend = id + nd1(1)*iintfp + uhf*nd1(2)*iintfp
      mxcor = maxcor - iend + 1
      call zero(icore(id),nd1(1)+uhf*nd1(2))
      if (mxcor .le. 0) call insmem('ONERDM',iend,maxcor)
      call densoo(icore(id),icore(iend),mxcor,uhf,singles)

      do spin = 1,uhf+1
        itmp = 0
        iboff = id + nd1(1)*iintfp*(spin-1)
        call aces_list_touch(1,nd1(spin),spin,160,itmp)
        call putlst(icore(iboff),1,1,1,spin,160)
      end do

C VV Part of density
C All off-diagonal elements are zeroed out
      id = 1
      iend = id + nd2(1)*iintfp + uhf*nd2(2)*iintfp
      mxcor = maxcor - iend + 1
      call zero(icore(id),nd2(1)+uhf*nd2(2))
      if (mxcor .le. 0) call insmem('ONERDM',iend,maxcor)
      call densvv(icore(id),icore(iend),mxcor,uhf,singles)

      do spin = 1,uhf+1
        itmp = 0
        iboff = id + nd2(1)*iintfp*(spin-1)
        call aces_list_touch(1,nd2(spin),2+spin,160,itmp)
        call putlst(icore(iboff),1,1,1,2+spin,160)
      end do

      return
      end

      subroutine dddiff(a,b,n)
      implicit none
      integer ii,n
      double precision a(n),b(n)

      do ii = 1,n
        write(6,*) 'ii=',ii,' diff=',abs(a(ii)-b(ii))
      end do

      return
      end

