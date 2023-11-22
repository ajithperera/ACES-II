      subroutine oombpt_dodiis(fock,dens,errnew,bmat,oocycle,totlen,
     &                         scrlen,spin,uhf)
      implicit none
C Common Blocks
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer iflags(100)
      common/flags/iflags
      integer ooUnit,cyclimit,dstart
      common/oombpt/ooUnit,cyclimit,dstart
C Input variables
      integer oocycle,scrlen,totlen,spin,uhf
C Input / Local variables
      double precision fock(scrlen),dens(scrlen)
C Pre-allocated Local variables
      double precision errnew(scrlen),bmat(cyclimit+1,cyclimit+1)
C Local variables
      integer numcyc,cyc,iloc,iloc2
      double precision zilch,one,fact,lintol,half
      data zilch,one,half,lintol /0.0d0,1.0d0,0.5d0,1.0d-40/

      numcyc=min(oocycle-dstart+1,iflags(20))

      if (oocycle .le. dstart) then
          call oombptlist(fock,totlen,'f',spin,1,'r')
      else
        call getrec(20,'JOBARC','OOMPTRPP',(cyclimit+1)**2*iintfp,bmat)
        call blkcpy2(bmat,cyclimit+1,cyclimit+1,dens,numcyc+1,numcyc+1,
     &               1,1)
        if (uhf .eq. 1) call dscal((numcyc+1)**2,half,dens,1)
        call eig(dens,errnew,1,numcyc+1,0)
        do cyc=1,numcyc+1
          iloc=(cyc-1)*(numcyc+2)+1
          iloc2=numcyc+(cyc-1)*(numcyc+1)+1
          if (abs(dens(iloc)) .le. lintol) then
            dens(cyc)=zilch
          else
            dens(cyc)=-errnew(iloc2)/dens(iloc)
          endif
        end do
        call dgemv('n',numcyc+1,numcyc+1,one,errnew,numcyc+1,dens,1,
     &             zilch,dens(numcyc+2),1)
        call dcopy(numcyc+1,dens(numcyc+2),1,errnew,1)
C Generate new fock extrapolant
        call zero(fock,totlen)
        do cyc=1,numcyc
          fact=errnew(cyc)
          call oombptlist(dens,totlen,'f',spin,cyc,'r')
          call daxpy(totlen,fact,dens,1,fock,1)
        end do
      endif

      return
      end
