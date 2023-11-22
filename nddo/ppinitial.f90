subroutine ppinitial(coef,expa,xyz,ppout,lone,ltwo,imax)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::lone,ltwo
integer,intent(in),dimension(:,:)::imax
double precision,intent(inout)::ppout
dist=(xyz(1,1)-xyz(1,2))**2+(xyz(2,1)-xyz(2,2))**2+(xyz(3,1)-xyz(3,2))**2

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
        
         tsum=0.0D0
         do 51 it=1,imax(1,1)
            do 52 jt=1,imax(1,2)
               r=0.0D0
               t=0.0D0
               ss=0.0D0
               r=expa(it,1)
               t=expa(jt,2)

           pii=(r*xyz(icart,1)+t*xyz(icart,2))/(r+t)
            piiai=pii-xyz(icart,1)  
           pjj=(r*xyz(jcart,1)+t*xyz(jcart,2))/(r+t)
            pjjbj=pjj-xyz(jcart,2)
           ss=((pi/(r+t))**(1.5))*exp((((-r)*t)/(r+t))*dist)
        rps=(piiai)*ss
        rpp=pjjbj*rps
       if(icart.eq.jcart)then
        rpp=rpp+(1.0D0/(2.0D0*(r+t)))*ss
        end if
        tsum=tsum+rpp*coef(it,1)*coef(jt,2)
 52   continue
 51   continue
      ppout=tsum
         return
         end

