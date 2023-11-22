      SUBROUTINE T2PPT2SQ(NO,NU,TI,O2,T2,VOE,V,o4,vo2,T2N,oeh,oep)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,E,F
      logical print
      common/flags/iflags(100)
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      common/timeinfo/timein,timenow,timetot,timenew
      DIMENSION TI(1),O2(Nu,NU,No,NO),T2(Nu,No,No,Nu),VOE(Nu,NU,No,NO),
     *T2N(NO,No,NU,Nu),vo2(nu,nu,no,no),v(1),o4(1),oeh(no),oep(nu)
      DATA ZERO/0.0D+0/,ONE/1.0D+0/TWO/2.0D+0/,FOUR/4.0D+0/,
     *     EIGHT/8.0D+0/,HALF/0.5D+0/
      print=iflags(1).gt.12
      call ro2hpp(0,no,nu,ti,o2)
c111111111111111111111111111111111111111111111111111111
      call rintpp(1,nu,voe)
      CALL TRANSQ( voe,NU)
      call vt4(no,nu,ti,voe,v,o4)
      CALL TRANSQ(O2,NOU)
      CALL INSITU(NU,NO,NO,NU,TI,O2,13)
      CALL TRANSQ(t2n,NOU)
      CALL INSITU(NU,NO,NO,NU,TI,t2n,13)
      CALL MATMULSK(O2,o4,T2N,NO2U,NU,NU,0,0)
      call timer(1)
      if(print)write(6,6517)timenew
c      call zclock('dgr 17  ',1)
 6517 format('Diagram no.17 required',f15.2,'  seconds')
c
c222222222222222222222222222222222222222222222222222222222222222222
c
      CALL TRANSQ( voe,NU)
      call tranmd(o2,no,no,nu,nu,34)
      call matmulsk(o2,voe,t2,no2u,nu,nu,1,0)
      call tranmd(o2,no,no,nu,nu,34)
      call insitu(no,no,nu,nu,ti,t2,13)
      call tranmd(t2,nu,no,no,nu,23)
c
      call insitu(no,no,nu,nu,ti,t2n,13)
      call transq(t2n,nou)
      call insitu(no,nu,nu,no,ti,t2n,13)
      do 302 i=1,nu
         call rpakv(35,nu,i,ti,ti(nu2+1),v)
         call trant3(v,nu,5)
         call matmulsk(v,t2(1,1,1,i),t2n,nu2,no2,nu,0,0)
 302  continue
      call insitu(nu,nu,no,no,ti,t2n,13)
      call transq(t2n,nou)
      call insitu(nu,no,no,nu,ti,t2n,13)
      call timer(1)
      if(print)write(6,6518)timenew
c      call zclock('dgr 18  ',1)
 6518 format('Diagram no.18 required',f15.2,'  seconds')
c44444444444444444444444444444444444444444444444444444444444444444444444
c
      call insitu(no,no,nu,nu,ti,o2,13)
      call transq(o2,nou)
      call rintvo(1,no,nu,ti,t2)
      CALL TRANMD(t2,NO,NU,NU,NO,1234)
      call insitu(no,nu,nu,no,ti,t2,12)
      call symt21(o2,no,nu,nu,no,23)
      CALL VECMUL(O2,NO2U2,HALF)
      call matmulsk(o2,t2,voe,nou,nou,nou,1,0)
      CALL VECMUL(O2,NO2U2,two)
      call desm21(o2,no,nu,nu,no,23)
      call rintvo(2,no,nu,ti,t2)
      call vecmul(t2,no2u2,half)
      call tranmd(t2,no,nu,nu,no,1234)
      call insitu(no,nu,nu,no,ti,t2,12)
      call matmulsk(o2,t2,voe,nou,nou,nou,0,1)
      call tranmd(o2,no,nu,nu,no,23)
      CALL TRANSQ(voe,NOU)
      CALL INSITU(NU,NO,NO,NU,TI,voe,13)
      do 303 i=1,nu
         call rpakv(35,nu,i,ti,ti(nu2+1),v)
      call matmulsk(voe,v,t2n(1,1,1,i),no2,nu,nu2,0,0)
 303  continue
c
      call matmulsk(o2,t2,voe,nou,nou,nou,1,1)
      CALL TRANSQ(voe,NOU)
      CALL INSITU(NU,NO,NO,NU,TI,voe,13)
      do 304 i=1,nu
         call rpakv(35,nu,i,ti,ti(nu2+1),v)
      call trant3(v,nu,2)
      call matmulsk(voe,v,t2n(1,1,1,i),no2,nu,nu2,0,0)
 304  continue
      call timer(1)
      if(print)write(6,6520)timenew
c      call zclock('dgr 20  ',1)
 6520 format('Diagram no.20 required',f15.2,'  seconds')
c55555555555555555555555555555555555555555555555555555555555555555555
      call rintvo(4,no,nu,ti,t2)
      call transq(t2,nou)
      call symt21(t2,nu,no,no,nu,23)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,13)
      call matmulsk(t2,o2,ti,nu,nu,no2u,1,1)
      call tranmd(o2,no,no,nu,nu,12)
      call matmulsk(o2,ti,t2n,no2u,nu,nu,0,0)
      CALL INSITU(NO,NO,NU,NU,TI,T2N,13)
      CALL TRANSQ(T2N,NOU)
      call timer(1)
      if(print)write(6,6521)timenew
c      call zclock('dgr 21  ',1)
 6521 format('Diagram no.21 required',f15.2,'  seconds')
c 666666666666666666666666666666666666666666666666666666666666666666
      call insitu(no,no,nu,nu,ti,o2,13)
      call transq(o2,nou)
      call insitu(no,nu,nu,no,ti,o2,13)
 112       continue
      call tranmd(o2,nu,nu,no,no,12)
      goto 1200
 1100 continue
      call ro2hpp(0,no,nu,ti,o2)
      call insitu(no,nu,nu,no,ti,o2,13)
      call tranmd(o2,nu,nu,no,no,12)
c
      call ro2hpp(1,no,nu,ti,voe)
      call t2den(voe,ti,oeh,oep,no,nu)
      call vecmul(voe,no2u2,half)
      call insitu(no,nu,nu,no,ti,voe,13)
      call veccop(no2u2,vo2,voe)
      call symt21(vo2,nu,nu,no,no,12)
c
      call ro2hpp(0,no,nu,ti,t2)
      call transq(t2,nou)
c
      call insitu(no,nu,nu,no,ti,t2n,13)
      do 133 ia=1,nu
         call rpakv(35,nu,ia,ti,ti(nu2+1),v)
      call trant3(v,nu,1)
         call zeroma(o4,1,nu3)
         do 132 i=1,no
         do 132 j=1,no
            call matmulsk(o2(1,1,i,j),v,ti,nu,nu2,nu,1,0)
            call trant3(ti,nu,3)   
            call matmulsk(voe(1,1,i,j),ti,o4,nu,nu2,nu,0,0)
            call trant3(ti,nu,1)   
            call matmulsk(voe(1,1,j,i),ti,o4,nu,nu2,nu,0,0)
            call trant3(ti,nu,3)   
            call matmulsk(vo2(1,1,j,i),ti,o4,nu,nu2,nu,0,1)
 132     continue
         call trant3(o4,nu,3)
         call matmulsk(o4,t2(1,1,1,ia),t2n,nu2,no2,nu,0,0)
 133  continue
      CALL INSITU(Nu,Nu,No,No,TI,T2N,13)
c
 1200 continue
      call timer(1)
      if(print)write(6,6522)timenew
c      call zclock('dgr 22  ',1)
 6522 format('Diagram no.22 required',f15.2,'  seconds')
c
c 777777777777777777777777777777777777777777777777777
c
      call tranmd(o2,nu,nu,no,no,12)
      call insitu(nu,nu,no,no,ti,o2,13)
      call rintvo(4,no,nu,ti,t2)
      call tranmd(t2,no,nu,nu,no,23)
      call transq(t2,nou)
      call insitu(nu,no,no,nu,ti,t2,12)
      call vecmul(t2,no2u2,two)
      call matmulsk(t2,o2,voe,nou,nou,nou,1,0)
      call vecmul(t2,no2u2,half)
      call tranmd(t2,no,nu,no,nu,13)
      call matmulsk(t2,o2,voe,nou,nou,nou,0,1)
      call tranmd(t2,no,nu,no,nu,13)
      call tranmd(o2,no,nu,nu,no,23)
      call matmulsk(t2,o2,voe,nou,nou,nou,0,1)
      call tranmd(o2,no,nu,nu,no,23)
      call insitu(no,nu,nu,no,ti,voe,12)
      call matmulsk(o2,voe,t2n,nou,nou,nou,0,0)
      call tranmd(o2,no,nu,nu,no,23)
      call vecmul(o2,no2u2,half)
      call matmulsk(o2,voe,t2n,nou,nou,nou,0,1)
      call vecmul(o2,no2u2,two)
      call vecmul(o2,no2u2,half)
      call tranmd(t2,no,nu,no,nu,13)
      call matmulsk(t2,o2,voe,nou,nou,nou,1,0)
      call vecmul(o2,no2u2,two)
      call tranmd(o2,no,nu,nu,no,23)
      call insitu(no,nu,nu,no,ti,voe,12)
      call matmulsk(o2,voe,t2n,nou,nou,nou,0,0)
      call tranmd(o2,no,nu,nu,no,23)
      call tranmd(t2n,no,nu,nu,no,23)
      call matmulsk(o2,voe,t2n,nou,nou,nou,0,0)
      call tranmd(o2,no,nu,nu,no,23)
      call tranmd(t2n,no,nu,nu,no,23)
      call timer(1)
      if(print)write(6,6523)timenew
c      call zclock('dgr 23  ',1)
 6523 format('Diagram no.23 required',f15.2,'  seconds')
c333333333333333333333333333333333333333333333333333333333333333333333
      call rintvo(1,no,nu,ti,voe)
      call transq(voe,nou)
      call insitu(nu,no,no,nu,ti,voe,13)
      do 311 i=1,nu
         call rpakv(35,nu,i,ti,ti(nu2+1),v)
      call trant3(v,nu,1)
      call trant3(v,nu,3)
      CALL MATMULSK(voe,v,t2(1,1,1,i),No2,Nu,NU2,1,0)
 311  continue
      call insitu(no,no,nu,nu,ti,t2,13)
      call transq(t2,nou)
      call insitu(no,nu,nu,no,ti,t2,13)
      call insitu(nu,nu,no,no,ti,t2,23)
      CALL TRANSQ(t2,NOU)
      call symt21(o2,no,nu,nu,no,23)
      CALL MATMULSK(o2,t2,T2N,NOU,NOU,NOU,0,0)
      call desm21(o2,no,nu,nu,no,23)
c 2
      do 312 i=1,nu
         call rpakv(35,nu,i,ti,ti(nu2+1),v)
      call trant3(v,nu,1)
      CALL MATMULSK(voe,V,t2(1,1,1,i),No2,Nu,NU2,1,0)
 312  continue
      call rintvo(2,no,nu,ti,voe)
      call transq(voe,nou)
      call insitu(nu,no,no,nu,ti,voe,13)
c 3
      do 313 i=1,nu
         call rpakv(35,nu,i,ti,ti(nu2+1),v)
      call trant3(v,nu,1)
      call trant3(v,nu,3)
      CALL MATMULSK(voe,V,t2(1,1,1,i),No2,Nu,NU2,0,0)
 313  continue
      call vecmul(voe,no2u2,-two)
c 4
      do 314 i=1,nu
         call rpakv(35,nu,i,ti,ti(nu2+1),v)
      call trant3(v,nu,1)
      CALL MATMULSK(voe,V,t2(1,1,1,i),No2,Nu,NU2,0,0)
 314  continue
c
      call insitu(no,no,nu,nu,ti,t2,13)
      call transq(t2,nou)
      call insitu(no,nu,nu,no,ti,t2,13)
c
      call insitu(nu,nu,no,no,ti,t2,23)
      CALL TRANSQ(t2,NOU)
      CALL MATMULSK(o2,t2,T2N,NOU,NOU,NOU,0,1)
      CALL TRANMD(o2,NO,NU,NU,NO,23)
      CALL TRANMD(t2n,NO,NU,NU,NO,23)
      CALL MATMULSK(o2,t2,T2N,NOU,NOU,NOU,0,1)
      CALL TRANMD(t2n,NO,NU,NU,NO,23)
      call timer(1)
      if(print)write(6,6519)timenew
      call zclock('t2ppt2sq',1)
 6519 format('Diagram no.19 required',f15.2,'  seconds')
      RETURN
      END
