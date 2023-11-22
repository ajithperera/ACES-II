subroutine dsover(coor,lone,ltwo,dsout,NGAUSS,coef,expa,rau)
use constants
implicit double precision (a-h,o-z)
integer:: dspot,sspot
double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::dsout
      dist=rau*rau
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
         rsum=0.0D0
         do  it=1,NGAUSS
            do  jt=1,NGAUSS
               r=0.0D0
               t=0.0D0
               ss=0.0D0
               r=expa(it,dspot)
               t=expa(jt,sspot)
               zeta=r+t
         pjj=(r*coor(jcart,dspot)+t*coor(jcart,sspot))/(zeta)
               pjjaj=pjj-coor(jcart,dspot)
         pii=(r*coor(icart,dspot)+t*coor(icart,sspot))/(zeta)
               piiai=pii-coor(icart,dspot)
        ss=((pi/(zeta))**(1.5))*exp((((-r)*t)/(zeta))*dist)
       rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(1.0D0/(2.0D0*zeta))*ss
          end if
           rsum=rsum+rds*coef(it,dspot)*coef(jt,sspot)
end do
end do
      dsout=rsum
      return
      end subroutine dsover
          
