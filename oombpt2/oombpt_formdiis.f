
      SUBROUTINE oombpt_formdiis(fock,dens,errnew,bmat,nmo,oocycle,
     &                           totlen,scrlen,spin)
      implicit none
C Common Blocks
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer iflags(100)
      common/flags/iflags
      integer ooUnit,cyclimit,dstart
      common/oombpt/ooUnit,cyclimit,dstart
C Input variables
      integer oocycle,scrlen,totlen,nmo,spin
C Input / Local variables
      double precision fock(scrlen),dens(scrlen)
C Pre-allocated Local variables
      double precision errnew(scrlen),bmat(cyclimit+1,cyclimit+1)
C Local variables
      integer nmo2,numcyc,thiscyc,cyc,bsize
      double precision zilch,onem,btmp,fact
      data onem,zilch /-1.0d0,0.0d0/
C External functions
      double precision ddot
C - - - - - - - - - - - - - - - -- - - - - - - - --  - - --  --- - - - - - - - -
C Basis set parameters
      nmo2=nmo*nmo
      if (iflags(20).gt.cyclimit) then
         print *, '@OOMBPTDIIS: Assertion failed.'
         print *, '          Requested ',iflags(20),' iterates but ',
     &            'only ',cyclimit,' are allowed.'
         call errex
      endif
      numcyc=min(oocycle-dstart+1,iflags(20))
      thiscyc=1+mod(oocycle-dstart,iflags(20))

C DIIS cycle
      bsize = (cyclimit+1)**2
      if (oocycle .lt. dstart) then
        call oombptlist(fock,totlen,'f',spin,1,'w')
      else if (oocycle .eq. dstart) then
C Setup DIIS
C Insert this fock matrix into full list
        call oombptlist(fock,totlen,'f',spin,1,'w')
C Construct error vector in original basis
        call form_err(fock,dens,errnew,totlen,spin)
C Store newest error vector in DIIS list
        call oombptlist(errnew,totlen,'d',spin,1,'w')
        if (spin .eq. 1) then
          call zero(bmat,bsize)
        else
          call getrec(20,'JOBARC','OOMPTRPP',bsize*iintfp,bmat)
        endif
        btmp=ddot(totlen,errnew,1,errnew,1)
        bmat(1,1)=btmp+bmat(1,1)
        call putrec(20,'JOBARC','OOMPTRPP',bsize*iintfp,bmat)
      else
C Actually use DIIS extrapolation
C Insert fock matrix into full list
        call oombptlist(fock,totlen,'f',spin,thiscyc,'w')
C Construct error vector in original basis
        call form_err(fock,dens,errnew,totlen,spin)
C at this point fock & dens are scratch
C Insert this error vector into DIIS list
        call oombptlist(errnew,totlen,'d',spin,thiscyc,'w')
C Form new contribution to B-matrix
        call getrec(20,'JOBARC','OOMPTRPP',bsize*iintfp,bmat)
        fact = dble(spin-1)
        do cyc=1,thiscyc-1
          call oombptlist(dens,totlen,'d',spin,cyc,'r')
          btmp=ddot(totlen,errnew,1,dens,1)
          bmat(cyc,thiscyc)=btmp+bmat(cyc,thiscyc)*fact
          bmat(thiscyc,cyc)=btmp+bmat(thiscyc,cyc)*fact
        end do
        btmp=ddot(totlen,errnew,1,errnew,1)
        bmat(thiscyc,thiscyc)=btmp+bmat(thiscyc,thiscyc)*fact
C Continue forming B matrix
        do cyc=thiscyc+1,numcyc
          call oombptlist(dens,totlen,'d',spin,cyc,'r')
          btmp=ddot(totlen,errnew,1,dens,1)
          bmat(cyc,thiscyc)=btmp+bmat(cyc,thiscyc)*fact
          bmat(thiscyc,cyc)=btmp+bmat(thiscyc,cyc)*fact
        end do
C at this point errnew is scratch
        do cyc=1,numcyc
          bmat(numcyc+1,cyc)=onem+bmat(numcyc+1,cyc)*fact
          bmat(cyc,numcyc+1)=onem+bmat(cyc,numcyc+1)*fact
          bmat(numcyc+1,numcyc+1)=zilch
        end do
        call putrec(20,'JOBARC','OOMPTRPP',bsize*iintfp,bmat)
      endif

      return
      end

      subroutine show_rot(rotmat,totlen,spin)
      implicit none
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop_full(8,2),vrt_full(8,2)
      common/sym_full/pop_full,vrt_full
      integer totlen,spin
      double precision rotmat(totlen),val,tol,one
      integer irrep,mos,ioff,ii,jj
      data tol,one /1.0d-3,1.0d0/

      ioff=1
      do irrep = 1,nirrep
        write(6,*) 'irrep=',irrep
        mos = pop_full(irrep,spin) + vrt_full(irrep,spin)
        do ii = 1,mos
          do jj = 1,mos
            val = rotmat(ioff+(ii-1)*mos+jj-1)
            if ((abs(val) .gt. tol) .and. (abs(val) .lt. (one-tol)))
     &        write(6,100) ii,jj,val
          end do
        end do
        ioff = ioff + mos*mos
      end do

      return
  100 format('ROTMAT(',I2,',',I2,')=',E16.10)
      end
