subroutine dsinitial(coef,expa,xyz,dsout,lone,ltwo,imax)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::lone,ltwo
integer,intent(in),dimension(:,:)::imax
double precision,intent(inout)::dsout
integer:: dspot,sspot
dist=(xyz(1,1)-xyz(1,2))**2+(xyz(2,1)-xyz(2,2))**2+(xyz(3,1)-xyz(3,2))**2

      if((lone.gt.4).and.(lone.lt.11))then
         dspot=1
         sspot=2
         if((lone.eq.5).or.(lone.eq.8).or.(lone.eq.9))then
             icart=1
            else if((lone.eq.6).or.(lone.eq.10))then
            icart=2
            else if(lone.eq.7)then
               icart=3
          end if          
          if(lone.eq.5)then
             jcart=1
          else if((lone.eq.6).or.(lone.eq.8))then
             jcart=2
          else if((lone.eq.7).or.(lone.eq.9).or.(lone.eq.10))then
             jcart=3
          end if
          end if


          if((ltwo.gt.4).and.(ltwo.lt.11))then
         dspot=2
         sspot=1
         if((ltwo.eq.5).or.(ltwo.eq.8).or.(ltwo.eq.9))then
             icart=1
            else if((ltwo.eq.6).or.(ltwo.eq.10))then
            icart=2
            else if(ltwo.eq.7)then
               icart=3
          end if          
          if(ltwo.eq.5)then
             jcart=1
          else if((ltwo.eq.6).or.(ltwo.eq.8))then
             jcart=2
          else if((ltwo.eq.7).or.(ltwo.eq.9).or.(ltwo.eq.10))then
             jcart=3
          end if
         end if
         w=0.0D0
         
         do 51 it=1,imax(1,dspot)
            do 52 jt=1,imax(1,sspot)
               r=0.0D0
               t=0.0D0
               ss=0.0D0
               r=expa(it,dspot)
               t=expa(jt,sspot)
               zeta=r+t
         pjj=(r*xyz(jcart,dspot)+t*xyz(jcart,sspot))/(zeta)
               pjjaj=pjj-xyz(jcart,dspot)
         pii=(r*xyz(icart,dspot)+t*xyz(icart,sspot))/(zeta)
               piiai=pii-xyz(icart,dspot)
        ss=((pi/(zeta))**(1.5))*exp((((-r)*t)/(zeta))*dist)
       rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(1.0D0/(2.0D0*zeta))*ss
          end if
           w=w+rds*coef(it,dspot)*coef(jt,sspot)
 52   continue
 51   continue
      dsout=w
   
      end subroutine dsinitial
          
