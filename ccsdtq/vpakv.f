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
