subroutine ppover(coor,lone,ltwo,ppout,NGAUSS,coef,expa,rau)
use constants
implicit double precision (a-h,o-z)
double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::ppout
dist=rau*rau
w=0.0D0
        if(lone.eq.2)then
           icart=1
        else if(lone.eq.3)then
           icart=2
         else if(lone.eq.4)then
            icart=3
         end if
            if(ltwo.eq.2)then
           jcart=1
        else if(ltwo.eq.3)then
           jcart=2
         else if(ltwo.eq.4)then
            jcart=3
            end if
  xxxx=coor(icart,1)
  yyyy=coor(icart,2)
  zzzz=coor(jcart,1)
  aaaa=coor(jcart,2)
            do it=1,NGAUSS
!            do  jt=1,NGAUSS
               r=expa(it,1)
               t=expa(1,2)
               rpt=r+t
           rpt=r+t
           pii=(r*xxxx+t*yyyy)/(rpt)
            piiai=pii-xxxx  
           pjj=(r*zzzz+t*aaaa)/(rpt)
            pjjbj=pjj-aaaa
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
        rps=(piiai)*ss
        rpp=pjjbj*rps
       if(icart.eq.jcart)then
        rpp=rpp+(1.0D0/(2.0D0*(rpt)))*ss
        end if
        w=w+rpp*coef(it,1)*coef(1,2)

               r=expa(it,1)
               t=expa(2,2)
               rpt=r+t
           rpt=r+t
           pii=(r*xxxx+t*yyyy)/(rpt)
            piiai=pii-xxxx  
           pjj=(r*zzzz+t*aaaa)/(rpt)
            pjjbj=pjj-aaaa
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
        rps=(piiai)*ss
        rpp=pjjbj*rps
       if(icart.eq.jcart)then
        rpp=rpp+(1.0D0/(2.0D0*(rpt)))*ss
        end if
        w=w+rpp*coef(it,1)*coef(2,2)


               r=expa(it,1)
               t=expa(3,2)
               rpt=r+t
           rpt=r+t
           pii=(r*xxxx+t*yyyy)/(rpt)
            piiai=pii-xxxx  
           pjj=(r*zzzz+t*aaaa)/(rpt)
            pjjbj=pjj-aaaa
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
        rps=(piiai)*ss
        rpp=pjjbj*rps
       if(icart.eq.jcart)then
        rpp=rpp+(1.0D0/(2.0D0*(rpt)))*ss
        end if
        w=w+rpp*coef(it,1)*coef(3,2)

               r=expa(it,1)
               t=expa(4,2)
               rpt=r+t
           rpt=r+t
           pii=(r*xxxx+t*yyyy)/(rpt)
            piiai=pii-xxxx  
           pjj=(r*zzzz+t*aaaa)/(rpt)
            pjjbj=pjj-aaaa
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
        rps=(piiai)*ss
        rpp=pjjbj*rps
       if(icart.eq.jcart)then
        rpp=rpp+(1.0D0/(2.0D0*(rpt)))*ss
        end if
        w=w+rpp*coef(it,1)*coef(4,2)

               r=expa(it,1)
               t=expa(5,2)
               rpt=r+t
           rpt=r+t
           pii=(r*xxxx+t*yyyy)/(rpt)
            piiai=pii-xxxx  
           pjj=(r*zzzz+t*aaaa)/(rpt)
            pjjbj=pjj-aaaa
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
        rps=(piiai)*ss
        rpp=pjjbj*rps
       if(icart.eq.jcart)then
        rpp=rpp+(1.0D0/(2.0D0*(rpt)))*ss
        end if
        w=w+rpp*coef(it,1)*coef(5,2)

               r=expa(it,1)
               t=expa(6,2)
               rpt=r+t
           rpt=r+t
           pii=(r*xxxx+t*yyyy)/(rpt)
            piiai=pii-xxxx  
           pjj=(r*zzzz+t*aaaa)/(rpt)
            pjjbj=pjj-aaaa
       ss=pi/rpt
       ss=ss*ss*ss
       ss=sqrt(ss)
       t1=-r*t/rpt
       t1=t1*dist
       t1=exp(t1)
       ss=ss*t1
        rps=(piiai)*ss
        rpp=pjjbj*rps
       if(icart.eq.jcart)then
        rpp=rpp+(1.0D0/(2.0D0*(rpt)))*ss
        end if
        w=w+rpp*coef(it,1)*coef(6,2)

 end do
      ppout=w
end subroutine ppover

