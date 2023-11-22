
      subroutine oombptdamp(fock,dens,oldvec,tmpvec,totlen,spin,ncycle)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop_full(8,2),vrt_full(8,2)
      common/sym_full/pop_full,vrt_full
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer ncycle,totlen,spin
C Input / Output variables
      double precision fock(totlen),dens(totlen)
C Pre-allocated Local variables
      double precision oldvec(totlen),tmpvec(totlen)
C External functions
      double precision dsum
      external dsum
C Local variables
      integer irrep,mos,ii,ioff
      double precision c,s,close_tol,lambda,densdiff,zilch,one,halfm
      character*1 sp(2)
      data close_tol,zilch,one,halfm /1.d-6,0.0d0,1.0d0,-0.5d0/
      data sp /'A','B'/

      if (ncycle .eq. 1) then
        call putrec(20,'JOBARC','OOPTDNS'//sp(spin),totlen*iintfp,dens)
        call putrec(20,'JOBARC','OOPTFCK'//sp(spin),totlen*iintfp,fock)
        return
      endif
      s=zilch
      c=zilch
      call getrec(20,'JOBARC','OOPTDNS'//sp(spin),totlen*iintfp,oldvec)

C Form density difference in the dens vector
      call daxpy(totlen,-one,oldvec,1,dens,1)
C Find largest density difference
      densdiff=zilch
      do ii = 1,totlen
        if (abs(densdiff) .lt. abs(dens(ii))) densdiff = abs(dens(ii))
      end do
C Form fock difference in fock vector
      call getrec(20,'JOBARC','OOPTFCK'//sp(spin),totlen*iintfp,oldvec)
      call daxpy(totlen,-one,oldvec,1,fock,1)
      ioff = 1
      do irrep = 1,nirrep
        mos = pop_full(irrep,spin) + vrt_full(irrep,spin)
        if (mos .gt. 0) then
          call xgemm('n','n',mos,mos,mos,one,oldvec(ioff),mos,
     &               dens(ioff),mos,zilch,tmpvec,mos)
          s = s+dsum(mos,tmpvec,mos+1)
          call xgemm('n','n',mos,mos,mos,one,fock(ioff),mos,dens(ioff),
     &               mos,zilch,tmpvec,mos)
          c = c + dsum(mos,tmpvec,mos+1)
        endif
        ioff = ioff + mos*mos
      end do

      if (densdiff .gt. close_tol) then
        if (c .lt. (s*halfm)) then
          lambda = one
        else
          lambda = s*halfm/c
        endif
        call dscal(totlen,lambda,fock,1)
        call dscal(totlen,lambda,dens,1)
      endif

      call daxpy(totlen,one,oldvec,1,fock,1)
      call getrec(20,'JOBARC','OOPTDNS'//sp(spin),totlen*iintfp,oldvec)
      call daxpy(totlen,one,oldvec,1,dens,1)

      call putrec(20,'JOBARC','OOPTDNS'//sp(spin),totlen*iintfp,dens)
      call putrec(20,'JOBARC','OOPTFCK'//sp(spin),totlen*iintfp,fock)

      return
      end
