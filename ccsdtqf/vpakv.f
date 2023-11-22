      subroutine vpakv(nfile,nu,ia,ib,ifr,inr,v)
      implicit double precision(a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/buforv/buf(341),ibuf(341),iiiii
      common/countv/irec
      dimension v(1),ifr(1),inr(1)
      data tr/0.1d-13/
      call zeroma(buf,1,341)
      call izo(ibuf,341)
      nrc=nu2/1024+1
      if(irec.eq.2*nrc)then
         call izo(ifr,nu2)
         call izo(inr,nu2)
      else
         call rifr(nfile,ifr,inr)
      endif
      ikol=(ia-1)*nu+ib
      irec=irec+1
      irc=irec
      ifr(ikol)=irec
      call wwp(nfile,nu2,irc,v)
      inr(ikol)=irc
      irec=irec+irc-1
      call wifr(nfile,ifr,inr)
      return
      end
      subroutine wwp(nfile,n,irc,v)
      implicit double precision(a-h,o-z)
      dimension v(n)
      common/buforv/buf(341),ibuf(341),iiiii
      data tr/0.1d-13/
      in=0
      nr=0
      irec=irc
      do 10 i=1,n
         x=v(i)
         if (dabs(x).gt.tr)then
            in=in+1
            buf(in)=x
            ibuf(in)=i
         endif
         if (in.eq.341)then
            nr=nr+1
            write(nfile,rec=irec)buf,ibuf,in
            in=0
            irec=irec+1
         endif
 10   continue
      write(nfile,rec=irec)buf,ibuf,in
      irc=nr+1
      return
      end


      subroutine rpakv(nfile,nu,ia,ifr,inr,v)
      implicit double precision(a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/buforv/buf(341),ibuf(341),iiiii
      dimension v(1),ifr(1),inr(1)
      call zeroma(v,1,nu3)
      call rifr(nfile,ifr,inr)
      do 100 i=1,nu
         ik=nu2*(i-1)
         ikol=(i-1)*nu+ia
         ifst=ifr(ikol)
         nrec=inr(ikol)
         kk=0
         do 40 ire=1,nrec
            irec=ifst+ire-1
            read(nfile,rec=irec)buf,ibuf,in
            do 41 ii=1,in
               k=ibuf(ii)+ik
               v(k)=buf(ii)
 41         continue
 40      continue
 100  continue
      return
      end
      subroutine rifr(nfile,ifr,inr)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension ifr(1),inr(1)
      nrec=nu2/1024+1
      mm=mod(nu2,1024)
      iof=1
      do 20 i=1,nrec-1
         call  ird(nfile,i,1024,ifr(iof))
         iof=iof+1024
 20   continue
      call  ird(nfile,nrec,mm,ifr(iof))
      iof=1
      do 30 i=1,nrec-1
         call  ird(nfile,i+nrec,1024,inr(iof))
         iof=iof+1024
 30   continue
      call  ird(nfile,2*nrec,mm,inr(iof))
      return
      end
      subroutine wifr(nfile,ifr,inr)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension ifr(1),inr(1)
      nrec=nu2/1024+1
      mm=mod(nu2,1024)
      iof=1
      do 20 i=1,nrec-1
         call  iwrt(nfile,i,1024,ifr(iof))
         iof=iof+1024
 20   continue
      call  iwrt(nfile,nrec,mm,ifr(iof))
      iof=1
      do 30 i=1,nrec-1
         call  iwrt(nfile,i+nrec,1024,inr(iof))
         iof=iof+1024
 30   continue
      call  iwrt(nfile,2*nrec,mm,inr(iof))
      return
      end

      subroutine rpakvt(nfile,nu,ia,ifr,inr,v,vcd)
      implicit double precision(a-h,o-z)
      integer a,b
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/buforv/buf(341),ibuf(341),iiiii
      dimension v(1),ifr(1),inr(1),vcd(nu,nu,nu)
      call rifr(nfile,ifr,inr)
      do 100 i=1,nu
         ikol=(i-1)*nu+ia
         ifst=ifr(ikol)
         nrec=inr(ikol)
         call rdvv(nfile,nu2,ifst,v)
         do 10 a=1,nu
            do 10 b=1,a
               iab=(a-1)*nu+b
               vcd(a,b,i)=v(iab)
 10      continue
         ikol=(ia-1)*nu+i
         ifst=ifr(ikol)
         nrec=inr(ikol)
         call rdvv(nfile,nu2,ifst,v)
         do 11 a=1,nu
            do 11 b=1,a
               iab=(a-1)*nu+b
               vcd(b,a,i)=v(iab)
 11      continue
 100  continue
      return
      end
      subroutine symvt(nu,ifr,inr,v,vcd)
      implicit double precision (a-h,o-z)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      dimension ifr(1),inr(1),v(1),vcd(nu,nu,nu)
      common/countv/irec
      nrec=nu2/1024+1
      irec=2*nrec
      do 10 ia=1,nu
         call rpakvt(36,nu,ia,ifr,inr,v,vcd)
         call trant3(vcd,nu,2)
         do 11 ib=1,nu
            call vpakv(35,nu,ib,ia,ifr,inr,vcd(1,1,ib))
 11      continue
 10   continue
      return
      end
      subroutine rdvv(nfile,n,ifst,v)
      implicit double precision (a-h,o-z)
      common/buforv/buf(341),ibuf(341),iiiii
      dimension v(1)
      call zeroma(v,1,n)
      ire=1
 10   continue
         irec=ifst+ire-1
         read(nfile,rec=irec)buf,ibuf,in
         do 41 ii=1,in
            k=ibuf(ii)
            v(k)=buf(ii)
 41      continue
         ire=ire+1
         if(in.eq.341)goto 10
         ifst=irec+1
      return
      end

      SUBROUTINE rdsymo4pak(nu,t2,t2s,svec,ifr,inr,nlist,npak)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION svec(1),t2(nu,nu),t2s(1),ifr(1),inr(1)
      integer svec,dirprd,a,b,c,d,pop,vrt,spop,svrt
      logical smlv
      character*8 lb
      common/countv/irec
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),N(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      common/sym/pop(8,2),vrt(8,2),ntaa(6)
      common/flags/iflags(100)
      equivalence(icllvl,iflags(2))
      smlv=icllvl.eq.11.or.icllvl.eq.13.or.icllvl.eq.14.or.
     *icllvl.eq.15.or.icllvl.eq.19.or.
     *icllvl.eq.21.or.icllvl.eq.22
      lb='SVAVB2  '
      call getrec(20,'JOBARC',lb,nu2,svec)
      kkk=0
      nrec=nu2/1024+1
      irec=2*nrec
      DO 10 IRREP=1,NIRREP
         NSYDIS=IRPDPD(IRREP,ISYTYP(2,233))
         NSYDSZ=IRPDPD(IRREP,ISYTYP(1,233))
         if(smlv)nsydsz=nsydis
         do 9 inum=1,nsydis
            call zeroma(t2,1,nu2)
            call rdlistnew(nlist,irrep,inum,1,t2s,nsydsz)
            do 8 is=1,nsydsz
               kkk=kkk+1
               call pckindsk22(kkk,3,svec,svec,nu2,nu2,a,b,c,d,233)
               t2(d,c)=t2s(is)
 8          continue
            call vpakv(npak,nu,a,b,ifr,inr,t2)
 9       CONTINUE
 10   CONTINUE
      call wifr(npak,ifr,inr)
      RETURN
      END
      subroutine iwrt(npak,i,n,ia)
      implicit double precision(a-h,o-z)
      common/buforv/ibuf(1024)
      dimension ia(n)
      call iveccop(n,ibuf,ia)
      write(npak,rec=i)ibuf
      return
      end
      subroutine ird(npak,i,n,ia)
      implicit double precision(a-h,o-z)
      common/buforv/ibuf(1024)
      dimension ia(n)
      read(npak,rec=i)ibuf
      call iveccop(n,ia,ibuf)
      return
      end
