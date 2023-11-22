

      subroutine oombpt_step_wrap(icore,maxcor,uhf,listv,ncycle)
      implicit none
C Parameters
      integer cyclim
      parameter (cyclim=25)
      integer dstart_param
      parameter (dstart_param=2)
C Common blocks
      integer ooUnit,cyclimit,dstart
      common/oombpt/ooUnit,cyclimit,dstart
C Common blocks
      integer iflags(100)
      common/flags/iflags
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
      integer nstart,nirrep
      common/syminf/nstart,nirrep
C Input variables
      integer uhf,listv,maxcor
C Input/output variables
      integer ncycle
C Pre-allocated local variables
      integer icore(maxcor)
C Local variables
      integer totlen,irrep,i0,i1,i2,i3,iend,maxcyc,ndrop,scrlen,nmo,
     &        mxirr2
      logical cycerror

C Didn't converge
      ncycle=ncycle+1
      cyclimit = cyclim
      dstart = dstart_param
      call form_dropmo_length(nmo,totlen,mxirr2)

      i0 = 1
      scrlen = max(totlen,(cyclimit+1)**2)
      i1=i0+iintfp*scrlen
      i2=i1+iintfp*scrlen
      i3=i2+iintfp*max(scrlen,nmo*nmo)
      iend=i3+iintfp*scrlen
      if ((iend-i0) .gt. maxcor)
     &                       call insmem('oombpt_step',iend-i0,maxcor)
C Form Fock and Density Matrix and do damping and diis
      call oombpt_step(icore(i0),icore(i1),icore(i2),icore(i3),totlen,
     &                 scrlen,nmo,mxirr2,listv,ncycle,uhf)

      call getrec(20,'JOBARC','OOMBPTCY',1,maxcyc)
      cycerror = ncycle .gt. maxcyc
      if (ncycle .le. maxcyc) then
        if (ncycle.lt.10) then
          write(6,444) ncycle
        else if (ncycle.lt.100) then
          write(6,445) ncycle
        else if (ncycle.lt.1000) then
          write(6,446) ncycle
        else
          write(6,447) ncycle
        endif
      endif
      call putrec(20,'JOBARC','OOMBPTIT',1,ncycle)
      call putrec(20,'JOBARC','OOMBPTST',1,0)

      call getrec(20,'JOBARC','IFLAGS',100,iflags)
C Reset parameters for running correlated 1-particle method
C Zero scf iterations
      iflags(16)=0
C If ROHF or UHF turn on UHF flag
C ROHF will be broken by the correlation
      if (iflags(11) .gt. 0) iflags(11)=1
C Turn on nonhf flag
      iflags(38)=1
C Make orbitals semicanonical
      iflags(39)=1
C Read guess from mos
      iflags(45)=0
C Turn off QRHF flags
      iflags(34)=0
      iflags(77)=0
      call putrec(20,'JOBARC','IFLAGS',100,iflags)
      call getrec(20,'JOBARC','NUMDROPA',1,ndrop)
C Reset some vector lengths
      if (ndrop .gt. 0) then
        call getrec(20,'JOBARC','MODROPA0',ndrop,icore)
        call putrec(20,'JOBARC','MODROPA ',ndrop,icore)
        call getrec(20,'JOBARC','REORDEA0',nmo,icore)
        call putrec(20,'JOBARC','REORDERA',nmo,icore)
        if (uhf .gt. 0) then
          call getrec(20,'JOBARC','MODROPB0',ndrop,icore)
          call putrec(20,'JOBARC','MODROPB ',ndrop,icore)
          call getrec(20,'JOBARC','REORDEB0',nmo,icore)
          call putrec(20,'JOBARC','REORDERB',nmo,icore)
        endif
      else
         call getrec(20,'JOBARC','OCCUPYA ',nirrep,icore)
         call putrec(20,'JOBARC','OCCUPYA0',nirrep,icore)
         if (uhf .gt. 0) then
           call getrec(20,'JOBARC','OCCUPYB ',nirrep,icore)
           call putrec(20,'JOBARC','OCCUPYB0',nirrep,icore)
         endif
      endif

      if (cycerror) then
        write(6,*) '@OOMBPT2: The OO-MBPT(2) equations did not converge'
        call putrec(20,'JOBARC','OOMBPTST',1,2)
        call getrec(20,'JOBARC','OOMBITER',(maxcyc+1)*iintfp,icore(i0))
        i1 = i0 + (maxcyc+1)*iintfp
        call getrec(20,'JOBARC','OOMBITR2',(maxcyc+1)*iintfp,icore(i1))
        call print_cycle(icore(i0),icore(i1),ncycle)
      endif

      return
  444 FORMAT(T3,'This is cycle ',I1,' of the OO-MBPT(2) equations.')
  445 FORMAT(T3,'This is cycle ',I2,' of the OO-MBPT(2) equations.')
  446 FORMAT(T3,'This is cycle ',I3,' of the OO-MBPT(2) equations.')
  447 FORMAT(T3,'This is cycle ',I4,' of the OO-MBPT(2) equations.')
      end


      subroutine oombpt_step(fock,dens,scr1,scr2,totlen,scrlen,nmo,
     &                       mxirr2,listv,ncycle,uhf)
      implicit none
C Common blocks
      integer ooUnit,cyclimit,dstart
      common/oombpt/ooUnit,cyclimit,dstart
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer totlen,scrlen,nmo,listv,ncycle,uhf,mxirr2
C Pre-allocated local variables
      double precision fock(scrlen),dens(scrlen),scr1(scrlen),
     &                 scr2(scrlen)
C Local variables
      logical bExist,damp
      integer reclen,nmo2,stat,spin
      double precision zilch
      character*8 cscfevcs(2)
      character*8 cscfevcsf(2)
      character*1 sp(2)
      data sp /'A','B'/
      data cscfevcs   /'SCFEVCAS','SCFEVCBS'/
      data cscfevcsf /'SCFEVCA0','SCFEVCB0'/
      data zilch /0.0d0/

      damp = .false.
      nmo2=nmo*nmo
C Open OO-MBPT(2) storage file
      ooUnit=15
      inquire(iolength=reclen) zilch
      reclen=reclen*nmo2
      if (ncycle .eq. 1) then
        open(unit=ooUnit,file='OOMBPT',form='UNFORMATTED',
     &       access='DIRECT',recl=reclen,status='REPLACE',err=666,
     &       iostat=stat)
        write(unit=ooUnit,rec=1,err=666,iostat=stat) 0
      else
        inquire(file='OOMBPT',exist=bExist,err=666,iostat=stat)
        if (bExist) then
          open(unit=ooUnit,file='OOMBPT',form='UNFORMATTED',
     &         access='DIRECT',recl=reclen,status='OLD',err=666,
     &         iostat=stat)
        else
          write(6,*) 'OOMBPT FILE DOES NOT EXIST'
          call errex
        endif
      endif

      do spin = 1,uhf+1
C Form fock matrix
        call getlst(scr1,1,1,1,spin,listv)
        call expand_vai(scr1,scr2,spin)
        call formfock_d(scr2,fock,scr1,totlen,nmo,spin)
C Form density matrix
        call form_dens(dens,totlen,spin)
C Back transform fock matrix and density matrix
        if (ncycle .gt. 1) then
C First undo semi-canonical transformation
C for fock
          call backtrans_semi(scr2,scr1,totlen,spin)
          call transform_obj(fock,scr1,scr2,'t','r',totlen,spin)
          call transform_obj(scr2,scr1,fock,'n','l',totlen,spin)
C for density
          call transform_obj(dens,scr1,scr2,'t','r',totlen,spin)
          call transform_obj(scr2,scr1,dens,'n','l',totlen,spin)
C Back transform fock into original basis
          call getrec(20,'JOBARC','OOMPROT'//sp(spin),totlen*iintfp,
     &                scr1)
          call transform_obj(fock,scr1,scr2,'t','r',totlen,spin)
          call transform_obj(scr2,scr1,fock,'n','l',totlen,spin)
C For density
          call transform_obj(dens,scr1,scr2,'t','r',totlen,spin)
          call transform_obj(scr2,scr1,dens,'n','l',totlen,spin)
        endif
C Call damping
        if (damp)
     &    call oombptdamp(fock,dens,scr1,scr2,totlen,spin,ncycle)
C Call diis
        call oombpt_formdiis(fock,dens,scr1,scr2,nmo,ncycle,totlen,
     &                       scrlen,spin)
      end do
c
      do spin = 1,uhf+1

        call oombpt_dodiis(fock,dens,scr1,scr2,ncycle,totlen,scrlen,
     &                     spin,uhf)
C At this point, the density is irrelevant, so it can be used as a temporary variable
C Form rotation
        call form_u_fock(fock,dens,totlen,spin)
        if (ncycle .gt. 1) then
C Back transform rotation
          call getrec(20,'JOBARC','OOMPROT'//sp(spin),totlen*iintfp,
     &                scr2)
          call transform_obj(dens,scr2,scr1,'t','l',totlen,spin)
        endif
C Store new unitary rotation on JOBARC
        call putrec(20,'JOBARC','OOMPROT'//sp(spin),totlen*iintfp,dens)
        if (ncycle .gt. 1) then
          call backtrans_semi(dens,scr2,totlen,spin)
          call transform_obj(scr1,scr2,dens,'t','l',totlen,spin)
        endif

        call getrec(20,'JOBARC',cscfevcs(spin),nmo2*iintfp,scr1)
C Now actually do the OO-MBPT(2) step
        call transform_orbs(scr1,dens,scr2,nmo,nmo,totlen,mxirr2,spin)

        call putrec(20,'JOBARC',cscfevcsf(spin),nmo2*iintfp,scr1)

      end do

      if (uhf .eq. 0)
     &  call putrec(20,'JOBARC',cscfevcsf(2),nmo2*iintfp,scr1)

      close(unit=ooUnit)
      return
666   call aces_io_error('oomptdiis',ooUnit,Stat)
      end
