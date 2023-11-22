      subroutine formfd(icore,maxcor,iuhf)

      implicit double precision (a-h,o-z)

      integer dirprd,pop,vrt
      dimension icore(maxcor)
      dimension icmo(2)

      common/info/nocco(2),nvrto(2)
      common/machsp/iintln,ifltln,iintfp,ialone,ibitwd
      common/syminf/nstart,nirrep,irrepa(255,2),dirprd(8,8)
      common/sympop/irpdpd(8,22),isytyp(2,500),njunk(18)
      common/bassym/nbas(8),nbasis,nbassq,nbastt
      common/offset/ioffc(8)
      common/pert/ntpert,npert(8),ipert(8),ixpert,iypert,izpert,
     &            ixypert,ixzpert,iyzpert,itransx,itransy,itransz,
     &            nucind
      common/sym/pop(8,2),vrt(8,2),nt(2),nf1(2),nf2(2)
      common/lsym/nlenq(8),nlent(8)

      ianti=0

      icmo(1)=1
      icmo(2)=icmo(1)+nbassq*iuhf*iintfp
      iscr=icmo(2)+nbasis*nbasis*iintfp
      ieva=iscr+nbasis*nbasis*iintfp
      ievb=ieva+iintfp*iuhf*nbasis
      ifd=ievb+iintfp*nbasiS
      mxscr=nbasis*nbasis

      lenfd=0
      lenfd2=0
      lensd=0
      lensd2=0

      do 10 irrep=1,nirrep
       lenfd=max(lenfd,nlent(irrep))
       lenfd2=max(lenfd2,nlenq(irrep))
       lensd=max(lensd,irpdpd(irrep,21),irpdpd(irrep,22))
       lensd2=max(lensd2,irpdpd(irrep,21),irpdpd(irrep,22),
     &            irpdpd(irrep,19),irpdpd(irrep,20),
     &            irpdpd(irrep,9),irpdpd(irrep,10))
10    continue

      isd=ifd+iintfp*lenfd*(iuhf+1)
      isd3=isd+iintfp*lensd*(iuhf+1)
      iend=isd3+iintfp*lensd2
      if(iend.ge.maxcor) call errex
c
c get eigen vectors and symmetry pack them
c
      call getrec(20,'JOBARC','SCFEVCA0',nbasis*nbasis*iintfp,
     &            icore(iscr))
      call symc(icore(iscr),icore(icmo(1)),nbasis,nbas,.false.,1)
      if(iuhf.eq.1) then
          call getrec(20,'joBARC','SCFEVCB0',
     &                 nbasis*nbasis*iintfp,icore(iscr))
          call symc(icore(iscr),icore(icmo(2)),nbasis,nbas,
     &               .false.,2)
      endif

c
c fill ioffc
c i don't know why this offset for
c but i am keeping it for time being

      IOFFC(1)=1
      DO 1 IRREP=1,NIRREP-1
       IOFFC(IRREP+1)=IOFFC(IRREP)+NBAS(IRREP)*NBAS(IRREP)
1     CONTINUE 

c read orbital energies

      call getrec(20,'JOBARC','SCFEVLA0',iintfp*nbasis,
     &            icore(ieva))
      if(iuhf.eq.1) then
         call getrec(20,'JOBARC','SCFEVLB0',iintfp*nbasis,
     &             icore(ievb))
      endif
      
      ifdexp=iend
      ifdmo=ifdexp+iintfp*(iuhf+1)*lenfd2
      IEND2=IFDMO+IINTFP*LENFD2*(IUHF+1)
      MXCOR2=MAXCOR-IEND2
       
      do 100 irrepp=1,nirrep
        np=npert(irrepp)
        nsp=np
        if(np.eq.0) go to 100
        do 150 ip=1,np

c read in the fock matrix derivative integrals
c list 102 contians that info done in vdint

       call getlst(icore(ifd),ip,1,1,irrepp,102)
c       call kkk(10,icore(ifd))
c get the overlap matrix derivaive 
c list 170 contains that info
c done in formsd
c in mo basis ij 

       call getlst(icore(isd),ip,1,1,irrepp,270)

c do the same for uhf case

       if(iuhf.eq.1) then
          ifd2=ifd+iintfp*nlent(irrepp)
          isd2=isd+irpdpd(irrepp,21)*iintfp
          call getlst(icore(ifd2),ip,1,1,irrepp,103)
c          call kkk(10,icore(ifd2))
          call getlst(icore(isd2),ip,1,1,irrepp,271)
       endif

c expand the matric into full matrix

       ifdexp2=ifdexp+iintfp*nlenq(irrepp)*iuhf
       call matexp(irrepp,nbas,icore(ifd),icore(ifdexp),ianti)
c       call kkk(16,icore(ifdexp))
       if(iuhf.eq.1) call matexp(irrepp,nbas,icore(ifd2),
     &                           icore(ifdexp2),ianti)

       noccsq=irpdpd(irrepp,21)
       nvrtsq=irpdpd(irrepp,19)
       novsq=irpdpd(irrepp,9)
c  compute the occupied-occupied block

       ioff=ifdmo
      call formij(irrepp,icore(ifdexp),icore(icmo(1)),nbassq,
     &             icore(iscr),mxscr,iuhf,icore(ioff),noccsq,
     &             .true.)

c  add -1/2 s(i,j)^CHI (F(i,i) + F(j,j)) to f(i,j)^chi

      call sxeval1(irrepp,icore(ioff),icore(ieva),icore(isd),
     &               pop(1,1),ianti)

      if(iuhf.eq.1) then
      call sxeval1(irrepp,icore(ioff+iintfp*noccsq),
     &                icore(ievb),icore(isd2),pop(1,2),
     &                ianti)
      end if

c  add s(ij) contribution to fij(chi)

      call sdinoo(irrepp,icore(ioff),icore(isd),icore(iend2),
     &               mxcor2,iuhf)

      call updmoi(np,irpdpd(irrepp,21),irrepp,276,0,0)
      call putlst(icore(ioff),ip,1,1,irrepp,276)

      if(iuhf.eq.1) then
        call updmoi(np,irpdpd(irrepp,22),irrepp,277,0,0)
        call putlst(icore(ioff+iintfp*noccsq),ip,1,1,irrepp,277)
      endif  

       ioff=ioff+iintfp*(iuhf*noccsq+irpdpd(irrepp,22))

c convert fock derivative from ao to mo  
c fai 
   
       call formai(irrepp,icore(ifdexp),icore(icmo(1)),nbassq,
     &             icore(iscr),mxscr,iuhf,icore(ioff),novsq,
     &             .true.)
cn       call kkk(irpdpd(irrepp,9),icore(ioff))
c fai^chi - sum(j,k) S(k,j)^chi * <aj||ik> 
c done is sdinvo
c ioff contains fai^chi
c isd contains s(k,j)^chi 
c iend2 we need for setting up icore in sdinvo
c i am not sure if it is right
      call sdinvo(irrepp,icore(ioff),icore(isd),icore(iend2),
     &               mxcor2,iuhf)

c list 174 contains info about overlap derivative 
c in mo basis
C S(a,i) type info

        call getlst(icore(isd3),ip,1,1,irrepp,274)

c  fai^chi - sum(j,k) S(k,j)^chi * <aj||ik> -S(a,i)^chi*E(i)
c done is sxeval2

        call sxeval2(irrepp,icore(ioff),icore(ieva),icore(isd3),
     &               pop(1,1),vrt(1,1))

c info need to put into the record
c set up the list using upmoi
c  Xai^chi
c should be same as list 174 S(a,j)^chi
       
       call updmoi(np,irpdpd(irrepp,9),irrepp,280,0,0)         
       call putlst(icore(ioff),ip,1,1,irrepp,280)
c       call kkk(irpdpd(irrepp,9),icore(ioff))
c  do the same for uhf case

       if(iuhf.eq.1) then
         call getlst(icore(isd3),ip,1,1,irrepp,275)
         call sxeval2(irrepp,icore(ioff+iintfp*novsq),icore(ievb),
     &                icore(isd3),pop(1,2),vrt(1,2))
         call updmoi(np,irpdpd(irrepp,10),irrepp,281,0,0) 
         call putlst(icore(ioff+iintfp*novsq),ip,1,1,irrepp,281)
       endif
150    continue
100   continue
      return
      end
