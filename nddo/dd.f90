subroutine ddover(coor,lone,ltwo,ddout,NGAUSS,coef,expa,rau)
use constants
implicit double precision (a-h,o-z)
integer:: dispot,dkspot
double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::ddout
       dist=rau*rau
       dispot=1
       dkspot=2
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
         if((ltwo.eq.5).or.(ltwo.eq.8).or.(ltwo.eq.9))then
             kcart=1
            else if((ltwo.eq.6).or.(ltwo.eq.10))then
            kcart=2
            else if(ltwo.eq.7)then
               kcart=3
          end if          
          if(ltwo.eq.5)then
             lcart=1
          else if((ltwo.eq.6).or.(ltwo.eq.8))then
             lcart=2
          else if((ltwo.eq.7).or.(ltwo.eq.9).or.(ltwo.eq.10))then
             lcart=3
          end if
             w=0.0D0
         rsum=0.0D0
         do  it=1,NGAUSS
            do  jt=1,NGAUSS
               r=0.0D0
               t=0.0D0
               ss=0.0D0
               r=expa(it,dispot)
               t=expa(jt,dkspot)
               zeta=r+t
         pll=(r*coor(lcart,dispot)+t*coor(lcart,dkspot))/(zeta)
               pllbl=pll-coor(lcart,dkspot)
         pkk=(r*coor(kcart,dispot)+t*coor(kcart,dkspot))/(zeta)
               pkkbk=pkk-coor(kcart,dkspot)
         pjj=(r*coor(jcart,dispot)+t*coor(jcart,dkspot))/(zeta)
               pjjaj=pjj-coor(jcart,dispot)
         pii=(r*coor(icart,dispot)+t*coor(icart,dkspot))/(zeta)
               piiai=pii-coor(icart,dispot)
        ss=((pi/(zeta))**(1.5))*exp((((-r)*t)/(zeta))*dist)
       rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(1.0D0/(2.0D0*zeta))*ss
          end if
          rdp=pkkbk*rds
          if(icart.eq.kcart)then
             rdp=rdp+(1.0D0/(2.0D0*zeta))*pjjaj*ss
             end if
             if(jcart.eq.kcart)then
                rdp=rdp+(1.0D0/(2.0D0*zeta))*piiai*ss
             end if
        rdd=pllbl*rdp
        if(icart.eq.lcart)then
           rps=(pjjaj)*ss
        rpp=pkkbk*rps
       if(jcart.eq.kcart)then
        rpp=rpp+(1.0D0/(2.0D0*(zeta)))*ss
        end if
        rdd=rdd+(1.0D0/(2.0D0*zeta))*rpp
        end if
        if(jcart.eq.lcart)then
           rps=(piiai)*ss
        rpp=pkkbk*rps
       if(icart.eq.kcart)then
        rpp=rpp+(1.0D0/(2.0D0*(zeta)))*ss
        end if
        rdd=rdd+(1.0D0/(2.0D0*zeta))*rpp
        end if
        if(kcart.eq.lcart)then
            rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(1.0D0/(2.0D0*zeta))*ss
          end if
          rdd=rdd+(1.0D0/(2.0D0*zeta))*rds
          end if
          rsum=rsum+rdd*coef(it,dispot)*coef(jt,dkspot)
 end do
 end do
      ddout=rsum
      return
      end
          
       
