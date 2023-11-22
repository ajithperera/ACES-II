      subroutine dpover(coor,lone,ltwo,dpout,NGAUSS,coef,expa,rau)
use constants
implicit double precision (a-h,o-z)
integer:: dispot,pkspot
double precision,intent(in):: rau
double precision,dimension(:,:),intent(in)::coef,expa,coor
integer,intent(in)::lone,ltwo,ngauss
double precision,intent(inout)::dpout

      dist=rau*rau
      if((lone.gt.4).and.(lone.lt.11))then
         dispot=1
         pkspot=2
         end if


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
         rsum=0.0D0
         do  it=1,NGAUSS
            
  
               r=expa(it,dispot)
               t=expa(1,pkspot)
               zeta=r+t
               zetainv=one/zeta
         pkk=(r*coor(kcart,dispot)+t*coor(kcart,pkspot))*zetainv
               pkkbk=pkk-coor(kcart,pkspot)
         pjj=(r*coor(jcart,dispot)+t*coor(jcart,pkspot))*zetainv
               pjjaj=pjj-coor(jcart,dispot)
         pii=(r*coor(icart,dispot)+t*coor(icart,pkspot))*zetainv
               piiai=pii-coor(icart,dispot)




        ss=pi*zetainv
        ss=ss*ss*ss
        ss=dsqrt(ss)
        rxt=-r*t*zetainv*dist
         ss=ss*dexp(rxt)
       rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(half*zetainv)*ss
          end if
          rdp=pkkbk*rds
          if(icart.eq.kcart)then
             rdp=rdp+(half*zetainv)*pjjaj*ss
             end if
             if(jcart.eq.kcart)then
                rdp=rdp+(half*zetainv)*piiai*ss
             end if
             rsum=rsum+rdp*coef(it,dispot)*coef(1,pkspot)

           r=expa(it,dispot)
               t=expa(2,pkspot)
               zeta=r+t
               zetainv=one/zeta
         pkk=(r*coor(kcart,dispot)+t*coor(kcart,pkspot))*zetainv
               pkkbk=pkk-coor(kcart,pkspot)
         pjj=(r*coor(jcart,dispot)+t*coor(jcart,pkspot))*zetainv
               pjjaj=pjj-coor(jcart,dispot)
         pii=(r*coor(icart,dispot)+t*coor(icart,pkspot))*zetainv
               piiai=pii-coor(icart,dispot)




        ss=pi*zetainv
        ss=ss*ss*ss
        ss=dsqrt(ss)
        rxt=-r*t*zetainv*dist
         ss=ss*dexp(rxt)
       rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(half*zetainv)*ss
          end if
          rdp=pkkbk*rds
          if(icart.eq.kcart)then
             rdp=rdp+(half*zetainv)*pjjaj*ss
             end if
             if(jcart.eq.kcart)then
                rdp=rdp+(half*zetainv)*piiai*ss
             end if
             rsum=rsum+rdp*coef(it,dispot)*coef(2,pkspot)

r=expa(it,dispot)
               t=expa(3,pkspot)
               zeta=r+t
               zetainv=one/zeta
         pkk=(r*coor(kcart,dispot)+t*coor(kcart,pkspot))*zetainv
               pkkbk=pkk-coor(kcart,pkspot)
         pjj=(r*coor(jcart,dispot)+t*coor(jcart,pkspot))*zetainv
               pjjaj=pjj-coor(jcart,dispot)
         pii=(r*coor(icart,dispot)+t*coor(icart,pkspot))*zetainv
               piiai=pii-coor(icart,dispot)




        ss=pi*zetainv
        ss=ss*ss*ss
        ss=dsqrt(ss)
        rxt=-r*t*zetainv*dist
         ss=ss*dexp(rxt)
       rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(half*zetainv)*ss
          end if
          rdp=pkkbk*rds
          if(icart.eq.kcart)then
             rdp=rdp+(half*zetainv)*pjjaj*ss
             end if
             if(jcart.eq.kcart)then
                rdp=rdp+(half*zetainv)*piiai*ss
             end if
             rsum=rsum+rdp*coef(it,dispot)*coef(3,pkspot)

r=expa(it,dispot)
               t=expa(4,pkspot)
               zeta=r+t
               zetainv=one/zeta
         pkk=(r*coor(kcart,dispot)+t*coor(kcart,pkspot))*zetainv
               pkkbk=pkk-coor(kcart,pkspot)
         pjj=(r*coor(jcart,dispot)+t*coor(jcart,pkspot))*zetainv
               pjjaj=pjj-coor(jcart,dispot)
         pii=(r*coor(icart,dispot)+t*coor(icart,pkspot))*zetainv
               piiai=pii-coor(icart,dispot)




        ss=pi*zetainv
        ss=ss*ss*ss
        ss=dsqrt(ss)
        rxt=-r*t*zetainv*dist
         ss=ss*dexp(rxt)
       rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(half*zetainv)*ss
          end if
          rdp=pkkbk*rds
          if(icart.eq.kcart)then
             rdp=rdp+(half*zetainv)*pjjaj*ss
             end if
             if(jcart.eq.kcart)then
                rdp=rdp+(half*zetainv)*piiai*ss
             end if
             rsum=rsum+rdp*coef(it,dispot)*coef(4,pkspot)

r=expa(it,dispot)
               t=expa(5,pkspot)
               zeta=r+t
               zetainv=one/zeta
         pkk=(r*coor(kcart,dispot)+t*coor(kcart,pkspot))*zetainv
               pkkbk=pkk-coor(kcart,pkspot)
         pjj=(r*coor(jcart,dispot)+t*coor(jcart,pkspot))*zetainv
               pjjaj=pjj-coor(jcart,dispot)
         pii=(r*coor(icart,dispot)+t*coor(icart,pkspot))*zetainv
               piiai=pii-coor(icart,dispot)




        ss=pi*zetainv
        ss=ss*ss*ss
        ss=dsqrt(ss)
        rxt=-r*t*zetainv*dist
         ss=ss*dexp(rxt)
       rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(half*zetainv)*ss
          end if
          rdp=pkkbk*rds
          if(icart.eq.kcart)then
             rdp=rdp+(half*zetainv)*pjjaj*ss
             end if
             if(jcart.eq.kcart)then
                rdp=rdp+(half*zetainv)*piiai*ss
             end if
             rsum=rsum+rdp*coef(it,dispot)*coef(5,pkspot)


r=expa(it,dispot)
               t=expa(6,pkspot)
               zeta=r+t
               zetainv=one/zeta
         pkk=(r*coor(kcart,dispot)+t*coor(kcart,pkspot))*zetainv
               pkkbk=pkk-coor(kcart,pkspot)
         pjj=(r*coor(jcart,dispot)+t*coor(jcart,pkspot))*zetainv
               pjjaj=pjj-coor(jcart,dispot)
         pii=(r*coor(icart,dispot)+t*coor(icart,pkspot))*zetainv
               piiai=pii-coor(icart,dispot)




        ss=pi*zetainv
        ss=ss*ss*ss
        ss=dsqrt(ss)
        rxt=-r*t*zetainv*dist
         ss=ss*dexp(rxt)
       rps=(piiai)*ss
       rds=pjjaj*rps
       if(icart.eq.jcart)then
          rds=rds+(half*zetainv)*ss
          end if
          rdp=pkkbk*rds
          if(icart.eq.kcart)then
             rdp=rdp+(half*zetainv)*pjjaj*ss
             end if
             if(jcart.eq.kcart)then
                rdp=rdp+(half*zetainv)*piiai*ss
             end if
             rsum=rsum+rdp*coef(it,dispot)*coef(6,pkspot)
 
              

















 end do
      dpout=rsum
      return
      end subroutine dpover
                
