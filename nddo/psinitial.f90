subroutine psinitial(coef,expa,xyz,psout,lone,ltwo,imax)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::lone,ltwo
integer,intent(in),dimension(:,:)::imax
double precision,intent(inout)::psout
integer:: pspot,sspot,smax,pmax
dist=(xyz(1,1)-xyz(1,2))**2+(xyz(2,1)-xyz(2,2))**2+(xyz(3,1)-xyz(3,2))**2
         w=0.0D0
tsum=0.0D0
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
         do 51 it=1,imax(1,pspot)
            do 52 jt=1,imax(1,sspot)
               r=0.0D0
               t=0.0D0
               ss=0.0D0
               r=expa(it,pspot)
               t=expa(jt,sspot)
               Pii=(r*xyz(icart,pspot)+t*xyz(icart,sspot))/(r+t)
               piiai=Pii-xyz(icart,pspot)
        ss=((pi/(r+t))**(1.5))*exp((((-r)*t)/(r+t))*dist)
       rps=(piiai)*ss
        tsum=tsum+rps*coef(it,pspot)*coef(jt,sspot)
 52   continue
 51   continue
      psout=tsum
        end subroutine psinitial
  

