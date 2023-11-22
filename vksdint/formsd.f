
      subroutine formsd(icore,maxcor,iuhf)
      
      implicit double precision (a-h,o-z)
      integer dirprd,pop,vrt
      logical scf,nonhf,field,geom,third,magn,spin

      dimension icore(maxcor)
      dimension icmo(2)

      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      common/bassym/nbas(8),nbasis,nbassq,nbastt
      common/offset/ioffc(8)
      common/pert/ntpert,npert(8),ipert(8),ixpert,iypert,izpert,
     &            iyzpert,ixzpert,ixypert,itransx,itransy,itransz,
     &            nucind
      common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)
      common/lsym/nlenq(8),nlent(8)
      common/pert2/field,geom,third,magn,spin

       ianti=0

      icmo(1)=1
      icmo(2)=icmo(1)+nbassq*iuhf*iintfp
      iscr=icmo(2)+nbasis*nbasis*iintfp
      isd=iscr+nbasis*nbasis*iintfp
      lensd=0
      do 10 irrep=1,nirrep
        lensd=max(lensd,nlent(irrep))
10    continue
      iend=isd+iintfp*nbasis*nbasis
      mxscr=nbasis*nbasis
      call getrec(1,'JOBARC','SCFEVCA0',nbasis*nbasis*iintfp,
     &            icore(iscr))
      call symc(icore(iscr),icore(icmo(1)),nbasis,nbas,.false.,1)
      if(iuhf.eq.1) then
       call getrec(20,'JOBARC','SCFEVCB0',
     &             nbasis*nbasis*iintfp,icore(iscr))
       Call symc(icore(iscr),icore(icmo(2)),nbasis,nbas,.false.,2)
      endif
      ioffc(1)=1
      do 1 irrep=1,nirrep-1
       ioffc(irrep+1)=ioffc(irrep)+nbas(irrep)*nbas(irrep)
1     continue 

      isdexp=iend
      isdmo=isdexp+nbasis*nbasis*iintfp
      DO 100 irrepp=1,nirrep
        np=npert(irrepp)
        if(np.eq.0) go to 100
         Do 150 ip=1,np
          call getlst(icore(isd),ip,1,1,irrepp,101)
          call matexp(irrepp,nbas,icore(isd),
     &                            icore(isdexp),ianti)
          noccsq=irpdpd(irrepp,21)
          nvrtsq=irpdpd(irrepp,19)
          novsq=irpdpd(irrepp,9)
          ioff=isdmo
          call formij(irrepp,icore(isdexp),icore(icmo(1)),nbassq,
     &             icore(iscr),mxscr,iuhf,icore(ioff),noccsq,
     &             .false.)
          call updmoi(np,irpdpd(irrepp,21),irrepp,270,0,0)
          call putlst(icore(ioff),ip,1,1,irrepp,270)
c          call kkk(irpdpd(irrepp,21),icore(ioff))
          if(iuhf.eq.1) then
            ioff2=ioff+iintfp*noccsq
            call updmoi(np,irpdpd(irrepp,22),irrepp,271,0,0)
            call putlst(icore(ioff2),ip,1,1,irrepp,271)
c            call kkk(irpdpd(irrepp,22),icore(ioff))
          endif
          ioff=ioff+iintfp*(iuhf*nvrtsq+irpdpd(irrepp,20))
          call formai(irrepp,icore(isdexp),icore(icmo(1)),nbassq,
     &             icore(iscr),mxscr,iuhf,icore(ioff),novsq,
     &             .false.)
c          call kkk(9,icore(ioff))
          call updmoi(np,irpdpd(irrepp,9),irrepp,274,0,0)
          call putlst(icore(ioff),ip,1,1,irrepp,274)
c          call kkk(irpdpd(irrepp,9),icore(ioff))
          if(iuhf.eq.1) then
             ioff2=ioff+iintfp*novsq
             call updmoi(np,irpdpd(irrepp,10),irrepp,275,0,0)
             call putlst(icore(ioff2),ip,1,1,irrepp,275)
c             call kkk(irpdpd(irrepp,10),icore(ioff))
          endif
150    continue
100   continue
      return
      end
