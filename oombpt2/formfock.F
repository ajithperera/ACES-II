
      subroutine formfock_wrap(icore,maxcor,uhf,listv)
      implicit none
C Common blocks
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer maxcor,uhf,listv
C Pre-allocated local variables
      integer icore(maxcor)
C Local variables
      integer naobas(8),nao,drop,mxirr2,irrep,ivai,ieval,ievcs1,ievcs2,
     &        ifock,iend,nmo,lenvo
      logical cycerror

      call getrec(20,'JOBARC','NBASTOT',1,nmo)
      call getrec(20,'JOBARC','NBASTOT',1,nao)
      call getrec(20,'JOBARC','NUMDROPA',1,drop)
      if (drop.ne.0) then
         call getrec(20,'JOBARC','NUMBASI0',nirrep,naobas)
      else
         call getrec(20,'JOBARC','NUMBASIR',nirrep,naobas)
      endif
      mxirr2=0
      do irrep = 1,nirrep
        mxirr2=max(mxirr2,naobas(irrep)*naobas(irrep))
      end do
      lenvo = max(nt(1),uhf*nt(2))

      ivai = 1
      ieval = ivai + lenvo*iintfp
      ievcs1 = ieval + nmo*iintfp
      ievcs2 = ievcs1 + nao*nmo*iintfp
      ifock = ievcs2 + mxirr2*iintfp
      iend = ifock + mxirr2*iintfp
      if (iend .gt. maxcor) call insmem('FORMFOCK',iend,maxcor)
      call formfock(icore(ifock),icore(ivai),icore(ieval),icore(ievcs1),
     &              icore(ievcs2),mxirr2,lenvo,nmo,nao,uhf,drop,naobas,
     &              cycerror,listv)

      if (cycerror) call conv_set(icore,maxcor,.false.)

      return
      end


      subroutine formfock(fock,vai,evals,scfevc,vecs,mxirr2,lenvo,nmo,
     &                    nao,uhf,drop,naobas,cycerror,listv)
      implicit none
C Common blocks
      integer nstart,nirrep
      common/syminf/nstart,nirrep
      integer pop(8,2),vrt(8,2),nt(2),nd1(2),nd2(2)
      common/sym/pop,vrt,nt,nd1,nd2
      integer nocco(2),nvrto(2)
      common/info/nocco,nvrto
      integer iflags(100)
      common/flags/iflags
      integer iflags2(500)
      common/flags2/iflags2
      integer iintln,ifltln,iintfp
      common/machsp/iintln,ifltln,iintfp
C Input variables
      integer nmo,nao,mxirr2,uhf,lenvo,drop,naobas(8),listv
      double precision vai(lenvo)
C Pre-allocated Local variables
      double precision fock(mxirr2),evals(nmo),scfevc(nao*nmo),
     &                 vecs(mxirr2)
C Output variables
      logical cycerror
C Local variables
      integer basis2,spin,iscfc,iscfr,ievalo,ivai,irrep,occ,virt,mos,
     &        aos,ii,aa,ncycle,maxcycle,ievalv
      double precision one,zilch,alpha
      character*8 cscfevl0(2),cscfevcs(2)
      data one,zilch /1.0d0,0.0d0/
      data cscfevl0 /'SCFEVLA0','SCFEVLB0'/
      data cscfevcs /'SCFEVCAS','SCFEVCBS'/

C Get pre-existing fock operator
      basis2 = nao*nmo
      alpha =0.0d0
      do spin=1,uhf+1
        call getrec(20,'JOBARC',cscfevl0(spin),nmo*iintfp,evals)
        call getrec(20,'JOBARC',cscfevcs(spin),basis2*iintfp,scfevc)
        call getlst(vai,1,1,1,spin,listv)
        iscfc=1
        iscfr=1
        ievalo=1
        ievalv=nocco(spin)+1
        ivai = 1
        do irrep=1,nirrep
          occ=pop(irrep,spin)
          virt=vrt(irrep,spin)
          mos=occ+virt
          aos=naobas(irrep)
          if (min(occ,virt) .gt. 0) then
c            call printq2(vai(ivai),virt,occ,'VAIMAT  ')
            call zero(fock,mos*mos)
C Reconstruct the full fock matrix
C Diagonal OO/VV part
            call dcopy(occ,evals(ievalo),1,fock,mos+1)
            call dcopy(virt,evals(ievalv),1,fock(occ*mos+occ+1),mos+1)
            do ii = 1,virt
              evals(ievalv+ii-1) = one
            end do
            call daxpy(virt,alpha,evals(ievalv),1,fock(occ*mos+occ+1),
     &                 mos+1)

C Add OV part Fock + potential part
            call blkcpy(vai(ivai),virt,occ,fock,mos,mos,occ+1,1)
C symmetrize ov block
            do ii=1,occ
              do aa=occ+1,mos
                fock((aa-1)*mos+ii)=fock((ii-1)*mos+aa)
              end do
            end do
c            call printq2(fock,mos,mos,'FOCKMAT ')
C Now diagonalize the fock operator to generate new orbitals
            call eig(fock,vecs,1,mos,1)
c            call printq2(vecs,mos,mos,'VECSMAT ')
C Form new orbitals
            call xgemm('n','n',aos,mos,mos,one,
     &                 scfevc((iscfc-1)*nao+iscfr),nao,vecs,mos,zilch,
     &                 fock,aos)
            call blkcpy(fock,aos,mos,scfevc,nao,nmo,iscfr,iscfc)
          endif
          iscfr=iscfr+aos
          iscfc=iscfc+mos
          ievalo = ievalo+occ
          ievalv = ievalv+virt
          ivai = ivai + virt*occ
        end do
C Put new orbitals in place
c        call printq2(scfevc,nao,nmo,'SCFEVC  ')
        call putrec(20,'JOBARC',cscfevcs(spin),basis2*iintfp,scfevc)
        if (uhf.eq.0) call putrec(20,'JOBARC',cscfevcs(2),
     &                        basis2*iintfp,scfevc)
      end do
      if (drop .gt. 0) then
        call getrec(20,'JOBARC','MODROPA0',drop,evals)
        call putrec(20,'JOBARC','MODROPA ',drop,evals)
        call getrec(20,'JOBARC','REORDEA0',nmo,evals)
        call putrec(20,'JOBARC','REORDERA',nmo,evals)
        if (uhf .gt. 0) then
          call getrec(20,'JOBARC','MODROPB0',drop,evals)
          call putrec(20,'JOBARC','MODROPB ',drop,evals)
          call getrec(20,'JOBARC','REORDEB0',nmo,evals)
          call putrec(20,'JOBARC','REORDERB',nmo,evals)
        endif
      endif

C Didn't converge
      call getrec(-1,'JOBARC','OOMBPTIT',1,ncycle)
      ncycle=ncycle+1
      call getrec(20,'JOBARC','OOMBPTCY',1,maxcycle)
      cycerror = ncycle .gt. maxcycle
      if (ncycle .le. maxcycle) then
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
      if (iflags(11).gt.0) iflags(11)=1
C Turn on nonhf flag
      iflags(38)=1
C Make orbitals semicanonical
      iflags(39)=1
C Turn off QRHF flags
      iflags(34)=0
      iflags(77)=0
      call putrec(20,'JOBARC','IFLAGS',100,iflags)

      return
  444 FORMAT(T3,'This is cycle ',I1,' of the OO-MBPT(2) equations.')
  445 FORMAT(T3,'This is cycle ',I2,' of the OO-MBPT(2) equations.')
  446 FORMAT(T3,'This is cycle ',I3,' of the OO-MBPT(2) equations.')
  447 FORMAT(T3,'This is cycle ',I4,' of the OO-MBPT(2) equations.')

      return
      end
