      subroutine T2HPT2SQ(NO,NU,FPP,TI,O2,T2,VO2,V,voe,vo,T2N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      INTEGER A,B,C,E,F
      logical print
      common/flags/iflags(100)
      COMMON/NEWTAPE/NT2T4,NT4INT
      common/vrinpi/ifirst(512)
      COMMON/ECCCONT/CONECC(20)
      common/timeinfo/timein,timenow,timetot,timenew
      COMMON/NN/NO2,NO3,NO4,NU2,NU3,NU4,NOU,NO2U,NO3U,NOU2,NOU3,NO2U2
      DIMENSION TI(1),O2(1),T2(NO,No,NU,Nu),VO2(NO,No,NU,Nu),
     *T2N(NO,No,NU,Nu),FPP(1),V(1),voe(no,no,nu,nu),vo(nu,no,no,nu)
      DATA ZERO/0.0D+0/,ONE/1.0D+0/TWO/2.0D+0/,FOUR/4.0D+0/,
     *     EIGHT/8.0D+0/,HALF/0.5D+0/
      call timer(0)
      print=iflags(1).gt.12
      call transq(t2n,nou)
      call insitu(nu,no,no,nu,ti,t2n,13)
      call ro2hpp(0,no,nu,ti,o2)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,13)
      call tranmd(o2,no,no,nu,nu,12)
      call ro2hpp(1,no,nu,ti,voe)
      call transq(voe,nou)
      call insitu(nu,no,no,nu,ti,voe,13)
      call tranmd(voe,no,no,nu,nu,12)
      call ro2hpp(2,no,nu,ti,vo)
      call transq(vo,nou)
      call insitu(nu,no,no,nu,ti,vo,13)
      call tranmd(vo,no,no,nu,nu,12)
      call rinthh(1,no,v)
c111111111111111111111111111111111111111111111111111111111111111
      call vecmul(vo,no2u2,two)
      call vecsub(vo,voe,no2u2)
      CALL MATMUL(v,vo,fpp,1,NU2,NO2,1,1)
      call vecadd(vo,voe,no2u2)
      call vecmul(vo,no2u2,half)
      call tranmd(o2,no,no,nu,nu,12)
      CALL MATMUL(O2,FPP,T2N,NO2U,NU,NU,0,0)
      CALL INSITU(NO,NO,NU,NU,TI,T2N,13)
      CALL TRANSQ(T2N,NOU)
      call timer(1)
      if(print)write(6,6501)timenew
      call zclock('dgr   1  ',1)
 6501 format('Diagram no. 1 required',f15.2,'  seconds')
c222222222222222222222222222222222222222222222222222222222222222222
      call tranmd(voe,no,no,nu,nu,12)
      call vecmul(voe,no2u2,two)
      call vecsub(voe,vo,no2u2)
      call transq(v,no)
      CALL MATMUL(v,voe,t2,no,NOU2,NO,1,1)
      call vecadd(voe,vo,no2u2)
      call vecmul(voe,no2u2,half)
      CALL INSITU(NO,NO,NU,NU,TI,T2,13)
      call transq(t2,nou)
      CALL INSITU(no,no,nu,nu,ti,O2,23)
      call tranmd(o2,no,nu,no,nu,13)
      CALL MATMUL(O2,T2,T2N,NOU,NOU,NOU,0,0)
cccccccccccccccccccccccccccccccccccccccccc
      CALL MATMUL(v,voe,t2,no,NOU2,NO,1,1)
      CALL INSITU(NO,NO,NU,NU,TI,T2,13)
      call transq(t2,nou)
      call tranmd(o2,no,nu,no,nu,13)
      CALL MATMUL(O2,T2,T2N,NOU,NOU,NOU,0,1)
c----------------------------------------------------
      CALL MATMUL(v,vo,t2,no,NOU2,NO,1,1)
      CALL INSITU(NO,NO,NU,NU,TI,T2,13)
      call transq(t2,nou)
      call tranmd(t2n,no,nu,nu,no,23)
      CALL MATMUL(O2,T2,T2N,NOU,NOU,NOU,0,1)
      call tranmd(t2n,no,nu,nu,no,23)
      call timer(1)
      if(print)write(6,6502)timenew
      call zclock('dgr   2  ',1)
 6502 format('Diagram no. 2 required',f15.2,'  seconds')
c333333333333333333333333333333333333333333333333333333333333333333333
      call insitu(no,nu,no,nu,ti,o2,12)
      call transq(o2,nou)
      call rinthh(2,no2,v)
      call mtrans(v,no,2)
c-----------------------------------------------na diagramy 3 i 4-----
      call symt21(v,no,no,no,no,13)
      CALL MATMUL(V,vo,t2,NO2,NU2,NO2,1,0)
      call desm21(v,no,no,no,no,13)
      CALL MATMUL(V,voe,t2,NO2,NU2,NO2,0,1)
      CALL INSITU(NO,NO,NU,NU,TI,T2,23)
      call tranmd(t2n,no,nu,nu,no,23)
      CALL MATMUL(t2,o2,t2n,NOU,NOU,NOU,0,0)
      call tranmd(t2n,no,nu,nu,no,23)
      call tranmd(o2,no,nu,nu,no,23)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call vecmul(v,no4,two)
      call mtrans(v,no,8)
      CALL MATMUL(v,voe,vo2,no2,nu2,no2,1,0)
      call vecmul(v,no4,half)
      CALL INSITU(NO,NO,NU,NU,TI,vo2,23)
      call vecadd(t2,vo2,no2u2)
      CALL MATMUL(t2,o2,t2n,NOU,NOU,NOU,0,0)
      call vecmul(o2,no2u2,half)
      call tranmd(o2,no,nu,nu,no,23)
      CALL MATMUL(vo2,o2,t2n,NOU,NOU,NOU,0,1)
      call vecmul(o2,no2u2,two)
      call timer(1)
      if(print)write(6,6503)timenew
      call zclock('dgr   3  ',1)
 6503 format('Diagram no. 3 required',f15.2,'  seconds')
c4444444444444444444444444444444444444444444444444444444444444444
      call insitu(no,no,nu,nu,ti,vo,13)
      call transq(vo,nou)
      call insitu(no,no,nu,nu,ti,voe,13)
      call transq(voe,nou)
      call mtrans(v,no,12)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,13)
      call tranmd(o2,no,no,nu,nu,12)
      CALL MATMUL(V,o2,t2,NO2,NU2,NO2,1,0)
      CALL INSITU(NO,NO,NU,NU,TI,T2,23)
      CALL MATMUL(t2,voe,t2n,NOU,NOU,NOU,0,1)
      call tranmd(t2,no,nu,no,nu,13)
      call vecmul(voe,no2u2,two)
      call vecsub(voe,vo,no2u2)
      CALL MATMUL(t2,voe,t2n,NOU,NOU,NOU,0,0)
      call vecadd(voe,vo,no2u2)
      call vecmul(voe,no2u2,half)
      call tranmd(t2,no,nu,no,nu,13)
      call tranmd(t2n,no,nu,nu,no,23)
      CALL MATMUL(t2,vo,t2n,NOU,NOU,NOU,0,1)
      call tranmd(t2n,no,nu,nu,no,23)
      call timer(1)
      if(print)write(6,6504)timenew
      call zclock('dgr   4  ',1)
 6504 format('Diagram no. 4 required',f15.2,'  seconds')
c555555555555555555555555555555555555555555555555555555555555
      call rintvo(1,no,nu,ti,t2)
      call transq(t2,nou)
      call insitu(no,nu,nu,no,ti,vo,13)
      call tranmd(vo,nu,nu,no,no,12)
      call tranmd(t2,nu,no,no,nu,23)
      call matmul(vo,t2,fpp,nu,nu,no2u,1,0)
      call rintvo(2,no,nu,ti,vo2)
      call transq(vo2,nou)
      call tranmd(vo2,nu,no,no,nu,23)
      call vecmul(vo2,no2u2,half)
      call vecsub(t2,vo2,no2u2)
      call insitu(no,nu,nu,no,ti,voe,13)
      call tranmd(voe,nu,nu,no,no,12)
      call vecmul(voe,no2u2,two)
      call matmul(voe,t2,fpp,nu,nu,no2u,0,1)
      call vecmul(voe,no2u2,half)
      call vecmul(vo2,no2u2,four)
      call matmul(vo,vo2,fpp,nu,nu,no2u,0,1)
      call insitu(no,no,nu,nu,ti,o2,13)
      call transq(t2n,nou)
      call matmul(o2,fpp,t2n,no2u,nu,nu,0,0)
      call transq(t2n,nou)
      call timer(1)
      if(print)write(6,6505)timenew
      call zclock('dgr   5  ',1)
 6505 format('Diagram no. 5 required',f15.2,'  seconds')
c666666666666666666666666666666666666666666666666666666666666666666
      call rintvo(1,no,nu,ti,v)
      call  tranmd(v,no,nu,nu,no,14)
      call rintvo(2,no,nu,ti,t2)
      call vecmul(t2,no2u2,half)
      call  tranmd(t2,no,nu,nu,no,14)
      call transq(o2,nou)
      call insitu(no,nu,nu,no,ti,o2,12)
      call insitu(nu,nu,no,no,ti,vo,23)
      call insitu(nu,nu,no,no,ti,voe,23)
      call vecmul(voe,no2u2,half)
      call tranmd(voe,nu,no,nu,no,13)
      call vecsub(vo,voe,no2u2)
      call matmul(v,vo,vo2,nou,nou,nou,1,1)      
      call vecadd(vo,voe,no2u2)
      call vecmul(voe,no2u2,four)
      call vecsub(voe,vo,no2u2)
      call matmul(t2,voe,vo2,nou,nou,nou,0,1)      
      call vecadd(voe,vo,no2u2)
      call vecmul(voe,no2u2,half)
c
      call  tranmd(vo2,no,nu,nu,no,14)
      call vecmul(o2,no2u2,two)
      call matmul(vo2,o2,t2n,nou,nou,nou,0,0)      
      call vecmul(o2,no2u2,half)
      call matmul(v,vo,vo2,nou,nou,nou,1,1)      
      call vecmul(voe,no2u2,two)
      call matmul(t2,voe,vo2,nou,nou,nou,0,1)      
      call vecmul(voe,no2u2,half)
      call  tranmd(vo2,no,nu,nu,no,14)
      call  tranmd(o2,nu,no,nu,no,13)
      call matmul(vo2,o2,t2n,nou,nou,nou,0,1)      
c
      call matmul(v,voe,vo2,nou,nou,nou,1,1)      
      call vecmul(vo,no2u2,two)
      call matmul(t2,vo,vo2,nou,nou,nou,0,1)      
      call vecmul(vo,no2u2,half)
      call  tranmd(vo2,no,nu,nu,no,14)
      call tranmd(t2n,no,nu,nu,no,23)
      call matmul(vo2,o2,t2n,nou,nou,nou,0,1)      
      call transq(t2n,nou)
      call insitu(nu,no,no,nu,ti,t2n,13)
      call timer(1)
      if(print)write(6,6506)timenew
      call zclock('dgr   6  ',1)
 6506 format('Diagram no. 6 required',f15.2,'  seconds')
c77777777777777777777777777777777777777777777777777777777777777777777777777
      call rintvo(1,no,nu,ti,t2)
      call insitu(no,nu,nu,no,ti,t2,13)
      call insitu(nu,no,nu,no,ti,vo,12)
      call transq(vo,nou)
      call insitu(nu,no,no,nu,ti,vo,13)
      call veccop(no2u2,o2,vo)
      call vecmul(vo,no2u2,two)
      call insitu(nu,no,nu,no,ti,voe,12)
      call transq(voe,nou)
      call insitu(nu,no,no,nu,ti,voe,13)
      call vecsub(vo,voe,no2u2)
      call rintvo(2,no,nu,ti,vo2)
      call insitu(no,nu,nu,no,ti,vo2,13)
      call veccop(no2u2,fpp,t2)
      call insitu(nu,nu,no,no,ti,t2,13)
      call transq(t2,nou)
      call insitu(nu,no,no,nu,ti,t2,13)
c
      call insitu(nu,nu,no,no,ti,vo2,13)
      call transq(vo2,nou)
      call insitu(nu,no,no,nu,ti,vo2,13)
c
      call insitu(no,no,nu,nu,ti,o2,13)
      call transq(o2,nou)
      call insitu(no,nu,nu,no,ti,o2,13)
c
      call insitu(no,no,nu,nu,ti,vo,13)
      call transq(vo,nou)
      call insitu(no,nu,nu,no,ti,vo,13)
c
      call ro2hpp(0,no,nu,v,ti)
      call tranmd(ti,no,nu,nu,no,23)
      call transq(ti,nou)
      call insitu(nu,no,no,nu,v,ti,13)
c
      call tranmd(t2n,no,no,nu,nu,34)
      do 300 i=1,nu
      call matmul(o2,t2(1,1,1,i) ,v,nu2,nu,no2,1,0)     
      call matmul(vo,vo2(1,1,1,i),v,nu2,nu,no2,0,1)     
      call trant3(v,nu,3)
      call matmul(fpp,voe(1,1,1,i),v,nu2,nu,no2,0,0)     
      call trant3(v,nu,1)
      call matmul(ti,v,t2n(1,1,1,i),no2,nu,nu2,0,0)     
 300  continue
c
      call ro2hpp(0,no,nu,ti,o2)
      call tranmd(o2,no,nu,nu,no,23)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,13)
      call tranmd(t2n,no,no,nu,nu,34)
c
      call insitu(nu,nu,no,no,ti,vo,13)
      call transq(vo,nou)
      call insitu(nu,no,no,nu,ti,vo,13)
c
      call vecadd(vo,voe,no2u2)
      call vecmul(vo,no2u2,half)
      call insitu(no,no,nu,nu,ti,t2n,13)
      call transq(t2n,nou)
      call tranmd(t2n,no,nu,nu,no,23)
      call timer(1)
      if(print)write(6,6507)timenew
      call zclock('dgr   7  ',1)
 6507 format('Diagram no. 7 required',f15.2,'  seconds')
c 88888888888888888888888888888888888888888888888888888888888
      call rintvo(1,no,nu,ti,t2)
      call transq(t2,nou)
      call insitu(no,no,nu,nu,ti,o2,13)
      call  transq(o2,nou)
      call tranmd(o2,no,nu,nu,no,23)
      call symt21(o2,no,nu,nu,no,23)
      call matmul(o2,t2,vo2,nou,nou,nou,1,0)
      call desm21(o2,no,nu,nu,no,23)
      call rintvo(2,no,nu,ti,t2)
      call transq(t2,nou)
      call matmul(o2,t2,vo2,nou,nou,nou,0,1)
      call insitu(no,nu,no,nu,ti,vo2,12)
      call transq(vo2,nou)
      call tranmd(vo2,no,nu,nu,no,14)
      call insitu(no,no,nu,nu,ti,vo,13)
      call transq(vo,nou)
      call insitu(no,nu,nu,no,ti,vo,12)
      call tranmd(t2n,no,nu,nu,no,23)
      call matmul(vo2,vo,t2n,nou,nou,nou,0,1)
      call tranmd(t2n,no,nu,nu,no,23)
      call insitu(no,no,nu,nu,ti,voe,13)
      call transq(voe,nou)
      call insitu(no,nu,nu,no,ti,voe,12)
      call tranmd(voe,nu,no,nu,no,13)
      call matmul(vo2,voe,t2n,nou,nou,nou,0,1)
c
      call tranmd(o2,no,nu,nu,no,23)
      call matmul(o2,t2,vo2,nou,nou,nou,1,1)
      call insitu(no,nu,no,nu,ti,vo2,12)
      call transq(vo2,nou)
      call tranmd(vo2,no,nu,nu,no,14)
      call vecmul(voe,no2u2,two)
      call vecsub(vo,voe,no2u2)
      call matmul(vo2,vo,t2n,nou,nou,nou,0,1)
      call vecadd(vo,voe,no2u2)
      call vecmul(voe,no2u2,half)
      call timer(1)
      if(print)write(6,6508)timenew
      call zclock('dgr   8  ',1)
 6508 format('Diagram no. 8 required',f15.2,'  seconds')
c9999999999999999999999999999999999999999999999999999999999999999999999999
      call rintvo(8,no,nu,ti,t2)
      call vecmul(t2,no2u2,-two)
      call rintvo(-7,no,nu,ti,t2)
      call rintvo(-6,no,nu,ti,t2)
      call vecmul(t2,no2u2,-one)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call insitu(no,nu,nu,no,ti,o2,12)
      call tranmd(o2,nu,no,nu,no,13)
      call symt21(o2,nu,no,nu,no,13)
      call tranmd(t2,no,nu,nu,no,23)
      call matmul(t2,o2,vo2,nou,nou,nou,1,0)
      call desm21(o2,nu,no,nu,no,13)
      call insitu(nu,no,nu,no,ti,o2,12)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,12)
      call matmul(o2,vo2,t2n,nou,nou,nou,0,0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call insitu(no,nu,no,nu,ti,o2,12)
      call transq(o2,nou)
      call insitu(no,nu,nu,no,ti,o2,12)
      call matmul(t2,o2,vo2,nou,nou,nou,1,0)
      call rintvo(8,no,nu,ti,t2)
      call tranmd(t2,no,nu,nu,no,23)
      call tranmd(o2,nu,no,nu,no,13)
      call matmul(t2,o2,vo2,nou,nou,nou,0,1)
      call tranmd(o2,nu,no,nu,no,13)
      call insitu(nu,no,nu,no,ti,o2,12)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,12)
      call tranmd(o2,no,nu,no,nu,13)
      call matmul(o2,vo2,t2n,nou,nou,nou,0,1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call rintvo(6,no,nu,ti,t2)
      call rintvo(-7,no,nu,ti,t2)
      call insitu(no,nu,no,nu,ti,o2,12)
      call transq(o2,nou)
      call insitu(no,nu,nu,no,ti,o2,12)
      call tranmd(t2,no,nu,nu,no,23)
      call matmul(t2,o2,vo2,nou,nou,nou,1,0)
      call insitu(nu,no,nu,no,ti,o2,12)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,12)
      call tranmd(t2n,no,nu,nu,no,23)
      call matmul(o2,vo2,t2n,nou,nou,nou,0,0)
      call tranmd(t2n,no,nu,nu,no,23)
      call timer(1)
      if(print)write(6,6509)timenew
      call zclock('dgr   9  ',1)
 6509 format('Diagram no. 9 required',f15.2,'  seconds')
c11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11  11 11 11  11 11 
      call rintvo(1,no,nu,ti,t2)
      call  tranmd(t2,no,nu,nu,no,14)
      call insitu(nu,no,nu,no,ti,vo,12)
      call transq(vo,nou)
      call insitu(nu,no,nu,no,ti,vo,12)
      call insitu(nu,no,nu,no,ti,voe,12)
      call transq(voe,nou)
      call insitu(nu,no,no,nu,ti,voe,12)
      call vecmul(voe,no2u2,half)
      call tranmd(voe,no,nu,no,nu,13)
      call vecsub(vo,voe,no2u2)
      call matmul(vo,t2,vo2,nou,nou,nou,1,1)      
      call vecadd(vo,voe,no2u2)
      call rintvo(2,no,nu,ti,t2)
      call vecmul(t2,no2u2,half)
      call  tranmd(t2,no,nu,nu,no,14)
      call vecmul(voe,no2u2,four)
      call vecsub(voe,vo,no2u2)
      call matmul(voe,t2,vo2,nou,nou,nou,0,1)      
      call vecadd(voe,vo,no2u2)
      call vecmul(voe,no2u2,half)
      call  tranmd(vo2,no,nu,nu,no,14)
      call insitu(no,nu,no,nu,ti,o2,12)
      call transq(o2,nou)
      call insitu(no,nu,nu,no,ti,o2,12)
      call tranmd(o2,nu,no,nu,no,13)
      call vecmul(vo2,no2u2,two)
      call matmul(vo2,o2,t2n,nou,nou,nou,0,0)      
      call rintvo(1,no,nu,ti,t2)
      call  tranmd(t2,no,nu,nu,no,14)
      call matmul(vo,t2,vo2,nou,nou,nou,1,1)      
      call rintvo(2,no,nu,ti,t2)
      call  tranmd(t2,no,nu,nu,no,14)
      call matmul(voe,t2,vo2,nou,nou,nou,0,1)      
      call  tranmd(vo2,no,nu,nu,no,14)
      call  tranmd(o2,nu,no,nu,no,13)
      call matmul(vo2,o2,t2n,nou,nou,nou,0,1)      
cccccccccccccccccccccccccccccccccccccccccccccc
      call rintvo(1,no,nu,ti,t2)
      call  tranmd(t2,no,nu,nu,no,14)
      call matmul(voe,t2,vo2,nou,nou,nou,1,1)      
c
      call rintvo(2,no,nu,ti,t2)
      call  tranmd(t2,no,nu,nu,no,14)
      call matmul(VO,t2,vo2,nou,nou,nou,0,1)      
      call  tranmd(vo2,no,nu,nu,no,14)
      call tranmd(t2n,no,nu,nu,no,23)
      call matmul(vo2,o2,t2n,nou,nou,nou,0,1)      
      call timer(1)
      if(print)write(6,6511)timenew
      call zclock('dgr  11  ',1)
 6511 format('Diagram no.11 required',f15.2,'  seconds')
c14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 14 
      call rintvo(1,no,nu,ti,t2)
      call insitu(no,nu,nu,no,ti,t2,13)
      call insitu(no,nu,no,nu,ti,vo,23)
      call insitu(no,nu,no,nu,ti,voe,23)
      call tranmd(voe,no,no,nu,nu,12)
      call matmul(vo,t2,v,no2,no2,nu2,1,0)     
      call mtrans(v,no,8)
      call matmul(voe,t2,v,no2,no2,nu2,0,0)     
      call mtrans(v,no,8)
      call rintvo(2,no,nu,ti,t2)
      call insitu(no,nu,nu,no,ti,t2,13)
      call vecmul(vo,no2u2,two)
      call vecsub(vo,voe,no2u2)
      call matmul(vo,t2,v,no2,no2,nu2,0,1)     
      call vecadd(vo,voe,no2u2)
      call vecmul(vo,no2u2,half)
      call mtrans(v,no,7)
      call insitu(nu,no,nu,no,ti,o2,12)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,13)
      call transq(t2n,nou)
      call insitu(nu,no,no,nu,ti,t2n,13)
      call matmul(v,o2,t2n,no2,nu2,no2,0,0)     
      call insitu(no,no,nu,nu,ti,t2n,13)
      call transq(t2n,nou)
      call tranmd(t2n,no,nu,nu,no,23)
      call timer(1)
      if(print)write(6,6514)timenew
      call zclock('dgr  14  ',1)
 6514 format('Diagram no.14 required',f15.2,'  seconds')
c15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 15 
      call rintvo(1,no,nu,ti,t2)
      call rintvo(2,no,nu,ti,vo2)
      call vecmul(vo2,no2u2,half)
      call vecsub(t2,vo2,no2u2)
      call insitu(no,no,nu,nu,ti,voe,13)
      call transq(voe,nou)
      call insitu(no,nu,nu,no,ti,voe,13)
      call tranmd(voe,nu,nu,no,no,12)
      call matmul(t2,voe,fpp,no,no,nou2,1,1)
      call vecadd(t2,vo2,no2u2)
      call vecmul(vo2,no2u2,two)
      call insitu(no,no,nu,nu,ti,vo,13)
      call transq(vo,nou)
      call insitu(no,nu,nu,no,ti,vo,13)
      call matmul(vo2,vo,fpp,no,no,nou2,0,1)
      call vecmul(fpp,no2,two)
      call matmul(t2,vo,fpp,no,no,nou2,0,0)
      call transq(fpp,no)
      call insitu(no,no,nu,nu,ti,o2,13)
      call transq(o2,nou)
      call tranmd(o2,no,nu,nu,no,23)
      call matmul(o2,fpp,t2n,nou2,no,no,0,0)
      call timer(1)
      if(print)write(6,6515)timenew
      call zclock('dgr  15  ',1)
 6515 format('Diagram no.15 required',f15.2,'  seconds')
c 16 1616 1616 1616 1616 1616 1616 1616 1616 1616 1616 1616 1616 16
      call rintpp(1,nu,v)
      call transq(v,nu)
      call insitu(nu,nu,no,no,ti,vo,13)
      call transq(vo,nou)
      call insitu(nu,no,no,nu,ti,vo,13)
      call vecmul(vo,no2u2,two)
      call insitu(nu,nu,no,no,ti,voe,13)
      call transq(voe,nou)
      call insitu(nu,no,no,nu,ti,voe,13)
      call tranmd(voe,no,no,nu,nu,12)
      call vecsub(vo,voe,no2u2)
      call matmul(vo,v,fpp,no2,1,nu2,1,1)
      call vecadd(vo,voe,no2u2)
      call vecmul(vo,no2u2,half)
      call matmul(o2,fpp,t2n,nou2,no,no,0,1)
      call timer(1)
      if(print)write(6,6516)timenew
      call zclock('dgr  16  ',1)
 6516 format('Diagram no.16 required',f15.2,'  seconds')
c 12 12 12  12 12 12  12 12 12  12 12 12  12 12 12  12 12 12  12 12 12 
      call insitu(no,no,nu,nu,ti,vo,23)
      call insitu(no,no,nu,nu,ti,voe,23)
      call vecmul(voe,no2u2,two)
      call vecsub(voe,vo,no2u2)
      call matmul(voe,v,t2,no2u,nu,nu,1,0)
      call matmul(t2,o2,t2n,nou,nou,nou,0,0)
      call vecadd(voe,vo,no2u2)
      call vecmul(voe,no2u2,half)
      call matmul(voe,v,t2,no2u,nu,nu,1,0)
      call tranmd(o2,no,nu,nu,no,23)
      call matmul(t2,o2,t2n,nou,nou,nou,0,1)
      call matmul(vo,v,t2,no2u,nu,nu,1,0)
      call tranmd(t2n,no,nu,nu,no,23)
      call matmul(t2,o2,t2n,nou,nou,nou,0,1)
      call tranmd(t2n,no,nu,nu,no,23)
      call timer(1)
      if(print)write(6,6512)timenew
      call zclock('dgr  12  ',1)
 6512 format('Diagram no.12 required',f15.2,'  seconds')
c 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 13 
      call rintvo(6,no,nu,ti,t2)
      call rintvo(-7,no,nu,ti,t2)
      call insitu(no,nu,nu,no,ti,t2,13)
      call transq(o2,nou)
      call insitu(nu,no,no,nu,ti,o2,13)
      call matmul(o2,t2,v,no2,no2,nu2,1,0)
      call rintvo(8,no,nu,ti,t2)
      call insitu(no,nu,nu,no,ti,t2,13)
      call tranmd(o2,no,no,nu,nu,12)
      call matmul(o2,t2,v,no2,no2,nu2,0,0)
      call tranmd(o2,no,no,nu,nu,12)
      call transq(t2n,nou)
      call insitu(nu,no,no,nu,ti,t2n,13)
      call matmul(v,o2,t2n,no2,nu2,no2,0,0)
      call insitu(no,no,nu,nu,ti,t2n,13)
      call transq(t2n,nou)
      call timer(1)
      if(print)write(6,6513)timenew
      call zclock('dgr  13  ',1)
 6513 format('Diagram no.13 required',f15.2,'  seconds')
c 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10
      call insitu(no,nu,no,nu,ti,vo,12)
      call tranmd(vo,nu,no,no,nu,23)
      call zeroma(t2,1,no2u2)
      do 301 i=1,nu
         irc=ifirst(i)
         call rdvv(46,nu3,irc,v)
         call symt21(v,nu,nu,nu,1,13)
         call matmul(v,vo(1,1,1,i),t2,nu2,no2,nu,0,1)
 301  continue
      call insitu(nu,nu,no,no,ti,t2,13)
      call insitu(no,no,nu,nu,ti,o2,23)
csak16
      call matmul(o2,t2,t2n,nou,nou,nou,0,1)
      call tranmd(o2,no,nu,no,nu,13)
      call tranmd(t2n,no,nu,nu,no,23)
      call matmul(o2,t2,t2n,nou,nou,nou,0,1)
      call tranmd(o2,no,nu,no,nu,13)
      call tranmd(t2n,no,nu,nu,no,23)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call insitu(no,nu,no,nu,ti,voe,12)
      call zeroma(t2,1,no2u2)
      call hpvvoe(0,no,nu,v,voe,t2)
      call insitu(nu,nu,no,no,ti,t2,13)
      call matmul(o2,t2,t2n,nou,nou,nou,0,1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call zeroma(t2,1,no2u2)
      call hpvvoe(1,no,nu,v,voe,t2)
      call insitu(nu,nu,no,no,ti,t2,13)
      call tranmd(o2,no,nu,no,nu,13)
      call tranmd(t2n,no,nu,nu,no,23)
      call matmul(o2,t2,t2n,nou,nou,nou,0,0)
      call tranmd(t2n,no,nu,nu,no,23)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call zeroma(t2,1,no2u2)
      call hpvvoe(2,no,nu,v,voe,t2)
      call insitu(nu,nu,no,no,ti,t2,13)
      call matmul(o2,t2,t2n,nou,nou,nou,0,0)
999   continue
      call timer(1)
      if(print)write(6,6510)timenew
      call zclock('t2hpt2sq',1)
 6510 format('Diagram no.10 required',f15.2,'  seconds')
      RETURN
      END
