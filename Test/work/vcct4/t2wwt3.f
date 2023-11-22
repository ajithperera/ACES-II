      subroutine t2wwt3(NO,NU,T1,TI,V,voe,ve,T2N,oeh,oep)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER A,B,C,E,F
      logical print
      common/flags/iflags(100)
      COMMON/NEWTAPE/NT2T4,NT4INT
      COMMON/ECCCONT/CONECC(20)
      common/timeinfo/timein,timenow,timetot,timenew
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION t1(1),TI(1),ve(1),OEH(NO),oep(nu),
     *V(NU,NU,NU,No),voe(no,nu,nu,no),t2n(no,nu,nu,no)
      DATA ZERO/0.0D+0/,ONE/1.0D+0/TWO/2.0D+0/,FOUR/4.0D+0/,
     *     EIGHT/8.0D+0/,HALF/0.5D+0/
      print=iflags(1).gt.12
      call ro2hpp(1,no,nu,ti,voe)
      call vecmul(voe,no2u2,half)
      call t2den(voe,ti,oeh,oep,no,nu)
      call timer(1)
c1111111111111111111111111111111111111111111111111111111111111111111
      call zeroma(t1,1,nou)
      call t1wt3int(no,nu,ti,t1,voe)
cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      call rdvem4(0,no,nu,ti,v)
      call matmulsk(t1,v,t2n,no,nou2,nu,0,0)
      call timer(1)
      if(print)write(6,6531)timenew
 6531 format('Diagram no.1p required',f15.2,'  seconds')
chhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
      call rdvem4(1,nu,no,ti,v)
      call tranmd(v,no,no,no,nu,12)
      call trt1(no,nu,ti,t1)
      call transq(t2n,nou)
      call matmulsk(t1,v,t2n,nu,no2u,no,0,1)
      call transq(t2n,nou)
      call timer(1)
      if(print)write(6,6532)timenew
 6532 format('Diagram no.1h required',f15.2,'  seconds')
c55555555555555555555555555555555555555555555555555555555555555555555
      call rdvem4(0,no,nu,ti,v)
      call symt21(v,nu,nu,nu,no,23)
      call insitu(no,nu,nu,no,ti,voe,13)
      call matmulsk(v,voe,t1,nu,no,nou2,1,0)
      call insitu(nu,nu,no,no,ti,voe,13)
      call transq(voe,nou)
      call insitu(nu,no,no,nu,ti,voe,13)
      call trt1(nu,no,ti,t1)
      call rdvem4(1,nu,no,ti,v)
      call symt21(v,no,no,no,nu,23)
      call matmulsk(v,voe,t1,no,nu,no2u,0,1)
      call trt1(no,nu,ti,t1)
      call vecmul(t1,nou,half)
      call insitu(no,nu,nu,no,ti,t2n,13)
      call tranmd(t2n,nu,nu,no,no,12)
      call t2ft3int(no,nu,t2n,ti,t1)
      call tranmd(t2n,nu,nu,no,no,12)
      call insitu(nu,nu,no,no,ti,t2n,13)
      call timer(1)
      if(print)write(6,6533)timenew
 6533 format('Diagram no.5ph required',f15.2,'  seconds')
c 33333333333333333333333333333333333333333333333333333333333333333333
cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      call rdvem4(0,no,nu,ti,v)
      call tranmd(v,nu,nu,nu,no,13)
      call matmulsk(voe,v,ve,no2,nou,nu2,1,0)
      call tranmd(ve,no,no,nu,no,14)
      call tranmd(ve,no,no,nu,no,12)
      call timer(1)
      if(print)write(6,6534)timenew
 6534 format('Diagram no.3p required',f15.2,'  seconds')
c4444444444444444444444444444444444444444444444444444444444444444444444
chhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
      call rdvem4(1,nu,no,ti,v)
      call insitu(no,no,nu,nu,ti,voe,13)
      call transq(voe,nou)
      call symt21(voe,no,nu,nu,no,23)
      call matmulsk(v,voe,ve,no2,nou,nou,0,0)
      call desm21(voe,no,nu,nu,no,23)
      call tranmd(v,no,no,no,nu,23)
      call matmulsk(v,voe,ve,no2,nou,nou,0,1)
      call tranmd(voe,no,nu,nu,no,23)
      call tranmd(ve,no,no,nu,no,14)
      call matmulsk(v,voe,ve,no2,nou,nou,0,1)
      call insitu(no,no,nu,no,ti,ve,23)
      call tranmd(ve,no,nu,no,no,13)
      call t2wt3h(no,nu,t2n,ti,ve)
      call timer(1)
      if(print)write(6,6535)timenew
 6535 format('Diagram no.4h required',f15.2,'  seconds')
c 33333333333333333333333333333333333333333333333333333333333333333333
c hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
      call tranmd(voe,no,nu,nu,no,23)
      call zeroma(ve,1,no3u)
      call rdvem3(1,nu,no,ti,ve)
      call tranmd(ve,no,no,nu,no,14)
      call insitu(no,nu,nu,no,ti,voe,13)
      call matmulsk(voe,ve,v,nu2,nou,no2,1,0)
      call tranmd(v,nu,nu,nu,no,13)
      inrec=1
      do 1 i=1,no
         irc=inrec
         call wwp(46,nu3,irc,v(1,1,1,i))
         inrec=inrec+irc
 1    continue
      nrec=inrec-1
      call timer(1)
      if(print)write(6,6536)timenew
 6536 format('Diagram no.3h required',f15.2,'  seconds')
c444444444444444444444444444444444444444444444444444444444444444444444444
cppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      call rdvem4(0,no,nu,ti,v)
      call insitu(nu,nu,no,no,ti,voe,23)
      call tranmd(voe,nu,no,nu,no,13)
      nread=1
      nwrit=1
      do 80 i=1,no
      irc=nread
      call rdvv(46,nu3,irc,ti)
      nread=irc
      call symt21(voe(1,1,1,i),nu,no,nu,1,13)
      call matmulsk(v,voe(1,1,1,i),ti,nu2,nu,nou,0,0)
      call desm21(voe(1,1,1,i),nu,no,nu,1,13)
      call tranmd(v,nu,nu,nu,no,13)
      call matmulsk(v,voe(1,1,1,i),ti,nu2,nu,nou,0,1)
      call tranmd(voe(1,1,1,i),nu,no,nu,1,13)
      call trant3(ti,nu,1)
      call matmulsk(v,voe(1,1,1,i),ti,nu2,nu,nou,0,1)
      call trant3(ti,nu,1)
      call tranmd(v,nu,nu,nu,no,13)
      irc=nwrit+nrec
      call wwp(46,nu3,irc,ti)
      nwrit=nwrit+irc
 80   continue
      nread=1+nrec
      do 3 i=1,no
         irc=nread
         call rdvv(46,nu3,irc,v(1,1,1,i))
         nread=irc
 3    continue
      call insitu(no,nu,nu,no,ti,t2n,13)
      call t2wt3p(no,nu,t2n,ti,v)
      call insitu(nu,nu,no,no,ti,t2n,13)
      call timer(1)
      if(print)write(6,6537)timenew
 6537 format('Diagram no.4p required',f15.2,'  seconds')
c666666666666666666666666666666666666666666666666666666666666666666
cpppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      call tranmd(voe,nu,no,nu,no,13)
      call insitu(nu,no,nu,no,ti,voe,12)
      call zeroma(ve,1,no3u)
      call vmwt3(no,nu,ve,ti,voe)
      call zeroma(v,1,nou3)
      call rdvem3(0,no,nu,ti,v)
      call tranmd(v,nu,nu,no,nu,14)
      call insitu(no,nu,nu,no,ti,t2n,13)
      call matmulsk(v,ve,t2n,nu2,no2,nou,0,1)
      call insitu(nu,nu,no,no,ti,t2n,13)
      call timer(1)
      if(print)write(6,6538)timenew
 6538 format('Diagram no.6p required',f15.2,'  seconds')
c 222222222222222222222222222222222222222222222222222222222222222222
chhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
      call rdvem3(1,nu,no,ti,v)
      call vecmul(v,no3u,two)
      call tranmd(ve,no,nu,no,no,13)
      call matmulsk(ve,v,t2n,nou,nou,no2,0,1)
      call vecmul(v,no3u,half)
      call tranmd(v,no,no,nu,no,14)
      call matmulsk(ve,v,t2n,nou,nou,no2,0,0)
      call tranmd(v,no,no,nu,no,14)
      call tranmd(ve,no,nu,no,no,14)
      call matmulsk(ve,v,t2n,nou,nou,no2,0,0)
      call tranmd(v,no,no,nu,no,14)
      call tranmd(t2n,no,nu,nu,no,23)
      call matmulsk(ve,v,t2n,nou,nou,no2,0,0)
      call tranmd(t2n,no,nu,nu,no,23)
      call timer(1)
      if(print)write(6,6539)timenew
 6539 format('Diagram no.2h required',f15.2,'  seconds')
c66666666666666666666666666666666666666666666666666666666666666666
chhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
      call insitu(no,nu,nu,no,ti,voe,13)
      call zeroma(v,1,nou3)
      call vewt3(no,nu,v,ti,voe)
      call tranmd(v,nu,nu,nu,no,13)
      call rdvem3(1,nu,no,ti,ve)
      call insitu(no,no,nu,no,ti,ve,13)
      call tranmd(ve,nu,no,no,no,34)
      call insitu(no,nu,nu,no,ti,t2n,13)
      call matmulsk(v,ve,t2n,nu2,no2,nou,0,1)
      call insitu(nu,nu,no,no,ti,t2n,13)
      inrec=1
      do 4 i=1,no
         irc=inrec
         call wwp(46,nu3,irc,v(1,1,1,i))
         inrec=inrec+irc
 4    continue
      nrec=inrec-1
      call timer(1)
      if(print)write(6,6540)timenew
 6540 format('Diagram no.6h required',f15.2,'  seconds')
c22222222222222222222222222222222222222222222222222222222222222222222
cpppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
      call rdvem3(0,no,nu,ti,v)
      call insitu(nu,nu,no,nu,ti,v,13)
      call tranmd(v,no,nu,nu,nu,24)
c
      nread=1
      nwrit=1
      do 90 i=1,no
      irc=nread
      call rdvv(46,nu3,irc,ti)
      nread=irc
      call trant3(ti,nu,3)
      call vecmul(ti,nu3,two)
      call matmulsk(v,ti,t2n(1,1,1,i),nou,nu,nu2,0,1)
      call vecmul(ti,nu3,half)
      call tranmd(v,no,nu,nu,nu,23)
      call matmulsk(v,ti,t2n(1,1,1,i),nou,nu,nu2,0,0)
      call tranmd(v,no,nu,nu,nu,23)
      call trant3(ti,nu,1)
      call matmulsk(v,ti,t2n(1,1,1,i),nou,nu,nu2,0,0)
      call tranmd(v,no,nu,nu,nu,23)
      call tranmd(t2n(1,1,1,i),no,nu,nu,1,23)
      call matmulsk(v,ti,t2n(1,1,1,i),nou,nu,nu2,0,0)
      call tranmd(v,no,nu,nu,nu,23)
      call tranmd(t2n(1,1,1,i),no,nu,nu,1,23)
 90   continue
      call timer(1)
      if(print)write(6,6541)timenew
 6541 format('Diagram no.2p required',f15.2,'  seconds')
      call symetr(t2n,no,nu)
      call engy1(no,nu,v,t2n,voe,' t2wwt3 ')
      call zeroma(t2n,1,no2u2)
 9999 continue
      RETURN
      END
