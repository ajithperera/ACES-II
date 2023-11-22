subroutine dpinitial(coef,expa,xyz,dpout,lone,ltwo,imax)
use constants
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::coef,expa,xyz
integer,intent(in)::lone,ltwo
integer,intent(in),dimension(:,:)::imax
double precision,intent(inout)::dpout
integer:: dispot,pkspot
dist=(xyz(1,1)-xyz(1,2))**2+(xyz(2,1)-xyz(2,2))**2+(xyz(3,1)-xyz(3,2))**2

      if((lone.gt.4).and.(lone.lt.11))then
         dispot=1
         pkspot=2
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
          if(ltwo.eq.2)then
             kcart=1
             else if(ltwo.eq.3)then
                kcart=2
              else if(ltwo.eq.4)then
                 kcart=3
                 end if
          end if

          if((ltwo.gt.4).and.(ltwo.lt.11))then
         dispot=2
         pkspot=1
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
          if(lone.eq.2)then
             kcart=1
             else if(lone.eq.3)then
                kcart=2
              else if(lone.eq.4)then
                 kcart=3
                 end if
             end if
             w=0.0D0
         tsum=0.0D0
         do 51 it=1,imax(1,dispot)
            do 52 jt=1,imax(1,pkspot)
               r=0.0D0
               t=0.0D0
               ss=0.0D0
               r=expa(it,dispot)
               t=expa(jt,pkspot)
               zeta=r+t
         pkk=(r*xyz(kcart,dispot)+t*xyz(kcart,pkspot))/(zeta)
               pkkbk=pkk-xyz(kcart,pkspot)
         pjj=(r*xyz(jcart,dispot)+t*xyz(jcart,pkspot))/(zeta)
               pjjaj=pjj-xyz(jcart,dispot)
         pii=(r*xyz(icart,dispot)+t*xyz(icart,pkspot))/(zeta)
               piiai=pii-xyz(icart,dispot)
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
             tsum=tsum+rdp*coef(it,dispot)*coef(jt,pkspot)
 52   continue
 51   continue
      dpout=tsum
      return
      end
                
