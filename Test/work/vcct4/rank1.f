      subroutine rank1(no,nu,ti,o2,vo2,v,oeh,oep)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/timeinfo/timein,timenow,timetot,timenew
      DIMENSION TI(1),O2(1),v(1),vo2(NO,NO,NU,NU),oeh(no),oep(nu)
      common/newio/npp1,nhh,npp,ntvoe,nto3,no1,nu1
      common/vrinpi/ifirst(512)
      DATA ZERO/0.0D+0/,ONE/1.0D+0/TWO/2.0D+0/,FOUR/4.0D+0/,
     *     EIGHT/8.0D+0/,HALF/0.5D+0/
      print=iflags(1).gt.12
      nno2=nu*(nu+1)/2
      lnhh=8*no2
      lnpp=8*nu2
      lvoe=8*nou2
      call ro2hpp(1,no,nu,ti,vo2)
      call vecmul(vo2,no2u2,half)
      call t2den(vo2,ti,oeh,oep,no,nu)
      call ro2hpp(0,no,nu,ti,o2)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,13)
      call tranmd(o2,no,no,nu,nu,12)
c
c fhh fhh fhh fhh fhh fhh fhh fhh fhh fhh fhh fhh fhh fhh fhh fhh fhh fhh
      call symt21(o2,no,no,nu,nu,12)
      call tranmd(o2,no,no,nu,nu,12)
C o2:hhpp,vo2:hpph: 
      call matmul(o2,vo2,v,no,no,nou2,1,0)
      call vrinthh(1,no,v)
c fpp fpp fpp fpp fpp fpp fpp fpp fpp fpp fpp fpp fpp fpp  fpp fpp  fpp fpp 
      call transq(vo2,nou)
c:vo2:phhp,o2:hhpp
      call matmul(vo2,o2,v,nu,nu,no2u,1,1)
      call vrintpp(1,nu,v)
c now voe voe voe voe voe voe voe voe voe voe voe voe voe  voe voe  voe voe 
      call transq(vo2,nou)
      call insitu(no,no,nu,nu,ti,o2,13)
      call transq(o2,nou)
      call insitu(no,nu,nu,no,ti,o2,12)
c:vo2:hpph,o2:phph 
      call matmul(vo2,o2,v,nou,nou,nou,1,0)
c
      call insitu(nu,no,nu,no,ti,o2,12)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,13)
      call desm21(o2,no,no,nu,nu,12)
      call insitu(no,no,nu,nu,ti,o2,13)
      call transq(o2,nou)
      call insitu(no,nu,nu,no,ti,o2,12)
c
      call tranmd(vo2,no,nu,nu,no,23)
      call matmul(vo2,o2,v,nou,nou,nou,0,1)
      call tranmd(v,no,nu,nu,no,1234)
      call vrintvo(1,no,nu,ti,v)
c voe2 voe2 voe2 voe2 voe2 voe2 voe2 voe2 voe voe voe voe voe  voe voe  voe 
      call tranmd(o2,nu,no,nu,no,13)
c:vo2:hpph,o2:phph 
      call matmul(vo2,o2,v,nou,nou,nou,1,1)
      call tranmd(v,no,nu,nu,no,1234)
      call vrintvo(2,no,nu,ti,v)
c voe3 voe3 voe3  voe3 voe3 voe3  voe3 voe3 voe3  voe3 voe3 voe3  voe3 voe3
      call tranmd(vo2,no,nu,nu,no,23)
c:vo2:hpph,o2:phph 
      call matmul(vo2,o2,v,nou,nou,nou,1,0)
      call tranmd(v,no,nu,nu,no,1234)
      call vrintvo(3,no,nu,ti,v)
c  vhhhh  vhhhh  vhhhh  vhhhh  vhhhh  vhhhh  vhhhh  vhhhh  vhhhh  vhhhh
      call insitu(no,nu,nu,no,ti,vo2,13)
      call insitu(nu,no,nu,no,ti,o2,12)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,13)
c vo2:pphh,o2:hhpp
      call matmul(o2,vo2,v,no2,no2,nu2,1,0)
      call vrinthh(2,no2,v)
c vpppp vppp vpppp vppp vpppp vppp vpppp vppp vpppp vppp vpppp vppp vpppp
cvo2:pphh,o2:hhpp
      call insitu(nu,nu,no,no,ti,vo2,13)
      call insitu(no,no,nu,nu,ti,o2,13)
c vo2:no,nu,nu,no, o2:nu,no,no,nu
      call transq(vo2,nou)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,vo2,13)
      call insitu(no,nu,nu,no,ti,o2,13)
      call tranmd(o2,nu,nu,no,no,12)
      call tranmd(vo2,no,no,nu,nu,12)
      ifirst(1)=1
c o2:pphh, vo2:hhpp
      do 101 i=1,nu
      call matmul(o2,vo2(1,1,1,i),v,nu2,nu,no2,1,0)
      call trant3(v,nu,5)
      irc=ifirst(i)
      call wwp(46,nu3,irc,v)
      if(i.ne.nu)ifirst(i+1)=ifirst(i)+irc
 101  continue
c voe4 voe4 voe4 voe4 voe4 voe4 voe4 voe4 voe4 voe4 voe4 voe4 voe4 voe4
      call insitu(no,no,nu,nu,ti,vo2,13)
      call tranmd(vo2,nu,no,no,nu,14)
cv:pppp,vo2:pphh
      call zeroma(o2,1,no2u2)
      do 102 i=1,nu
         call rpakv(35,nu,i,ti,ti(nu2+1),v)
         call matmul(v,vo2(1,1,1,i),o2,nu2,no2,nu,0,0)
 102  continue
 398  format(5f14.10)
      call tranmd(vo2,nu,no,no,nu,14)
      call transq(vo2,nou)
      call insitu(no,nu,nu,no,ti,vo2,13)
      call insitu(nu,nu,no,no,ti,o2,13)
      call vrintvo(4,no,nu,ti,o2)
c voe5 voe5 voe5  voe5 voe5 voe5  voe5 voe5 voe5  voe5 voe5 voe5 voe5
      call rdov4(1,nu,no,ti,v)
      call matmul(vo2,v,o2,nu2,no2,no2,1,0)
      call insitu(nu,nu,no,no,ti,o2,13)
      call vrintvo(5,no,nu,ti,o2)
c voe6 voe6 voe6 voe6 voe6 voe6 voe6 voe6 voe6 voe6 voe6 voe6 voe6 voe6 
      call insitu(nu,nu,no,no,ti,vo2,13)
      call tranmd(vo2,no,nu,nu,no,23)
      call ro2hpp(1,no,nu,ti,o2)
      call insitu(no,nu,nu,no,ti,o2,12)
      call symt21(vo2,no,nu,nu,no,23)
      call matmul(vo2,o2,v,nou,nou,nou,1,0)
      call desm21(vo2,no,nu,nu,no,23)
      call vrintvo(6,no,nu,ti,v)
      call ro2hpp(2,no,nu,ti,o2)
      call insitu(no,nu,nu,no,ti,o2,12)
      call matmul(vo2,o2,v,nou,nou,nou,1,1)
      call vrintvo(7,no,nu,ti,v)
      call tranmd(vo2,no,nu,nu,no,23)
      call matmul(vo2,o2,v,nou,nou,nou,1,1)
      call vrintvo(8,no,nu,ti,v)
      call timer(1)
      if(print)write(6,104)timenew
 104  format('Rank1 routine required',f10.1,' seconds')
      END
