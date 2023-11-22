subroutine psover(coor,lone,ltwo,psout,NGAUSS,coef,expa,rau)
use constants
implicit double precision (a-h,o-z)

double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::psout

integer pspot,sspot

dist=rau*rau

if((lone.gt.1).and.(lone.lt.5))then
  pspot=1
  sspot=2
           if(lone.eq.2)then
           icart=1
           else if(lone.eq.3)then
           icart=2
           else if(lone.eq.4)then
           icart=3
           end if
endif
           
if((ltwo.gt.1).and.(ltwo.lt.5))then
pspot=2
sspot=1
             if(ltwo.eq.2)then
             icart=1
             else if(ltwo.eq.3)then
             icart=2
             else if(ltwo.eq.4)then
             icart=3
             end if
 end if


w=zero
xxxx=coor(icart,pspot)
yyyy=coor(icart,sspot)

do  it=1,NGAUSS
!do  jt=1,NGAUSS
       r=expa(it,pspot)
       t=expa(1,sspot)
       rpt=r+t
       Pii=(r*xxxx+t*yyyy)/(rpt)
       piiai=Pii-xxxx
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
       rps=(piiai)*ss
        w=w+rps*coef(it,pspot)*coef(1,sspot)

       r=expa(it,pspot)
       t=expa(2,sspot)
       rpt=r+t
       Pii=(r*xxxx+t*yyyy)/(rpt)
       piiai=Pii-xxxx
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
       rps=(piiai)*ss
        w=w+rps*coef(it,pspot)*coef(2,sspot)

       r=expa(it,pspot)
       t=expa(3,sspot)
       rpt=r+t
       Pii=(r*xxxx+t*yyyy)/(rpt)
       piiai=Pii-xxxx
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
       rps=(piiai)*ss
        w=w+rps*coef(it,pspot)*coef(3,sspot)

       r=expa(it,pspot)
       t=expa(4,sspot)
       rpt=r+t
       Pii=(r*xxxx+t*yyyy)/(rpt)
       piiai=Pii-xxxx
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
       rps=(piiai)*ss
        w=w+rps*coef(it,pspot)*coef(4,sspot)

       r=expa(it,pspot)
       t=expa(5,sspot)
       rpt=r+t
       Pii=(r*xxxx+t*yyyy)/(rpt)
       piiai=Pii-xxxx
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
       rps=(piiai)*ss
        w=w+rps*coef(it,pspot)*coef(5,sspot)

       r=expa(it,pspot)
       t=expa(6,sspot)
       rpt=r+t
       Pii=(r*xxxx+t*yyyy)/(rpt)
       piiai=Pii-xxxx
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
       rps=(piiai)*ss
        w=w+rps*coef(it,pspot)*coef(6,sspot)



end do
!end do
 
      psout=w
      return
      end
 

