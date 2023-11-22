subroutine local(twoelcl,ni,nj,d1a,d2a,d1b,d2b,p0a,p0b,p1a,p1b,p2a,p2b,rsq,r) 
use constants
implicit double precision (a-h,o-z)

 double precision,intent(in)::d1a,d2a,d1b,d2b,p0a,p0b,p1a,p1b,p2a,p2b,rsq,r
 integer, intent(in)::ni,nj
 double precision , dimension(22), intent(inout) ::twoelcl
!     p = P-pi and o = P-sigma
!     order of storage in twoelcl (stands for 2 electron local coordinates)
!     twoelcl(1)=(ss|ss)  twoelcl(2)=(ss|pp) twoelcl(3)=(ss|oo)
!     twoelcl(4)=(ss|so) twoelcl(5)=(pp|ss) twoelcl(6)=(oo|ss)
!     twoelcl(7)=(so|ss)xx twoelcl(8)=(so|pp)xx twoelcl(9)=(pp|so)
!     twoelcl(10)=(so|oo)xx twoelcl(11)=(oo|so) twoelcl(12)=(sp|sp)xx
!     13 = (so|so)xx 14 = (sp|po)xx 15 = (op|sp) 16 = (po|po)
!     17 = (p p | p p) 18 = (pp|p*p*) 19 = (pp|oo) 20 = (oo|pp)
!     21 = ( o o | o o) 22 = (pp*|pp*)


!     hydrogen-hydrogen
      if(ni==1.and.nj==1)then
       xx=p0a+p0b
!ccc         a00=(p0a+p0b)*(p0a+p0b)
        a00=(xx)*(xx)
         qq=one/dsqrt(rsq+a00)
         twoelcl(1)=qq

!     hydrogen - heavy atom+++++++++++++++++++++++++++++++++++++++++++=


      elseif((ni==1).and.(nj==4))then
       
!    there are four such integrals as follows
!     (s s| s s) = [q,q]
!     (s s | p p) = [q,q] + [q,Qpp]
!     (s s | o o) = [q,q] + [q,oo]
!     (s s | o s) = [q,uo]


!     first the (s s | s s) integral
      xx=p0a+p0b
!ccc      a00=(p0a+p0b)*(p0a+p0b)
      a00=(xx)*(xx)
      qq=one/dsqrt(rsq+a00)
      twoelcl(1)=qq
!******************************

!   now the (s s | p p ) integral
      xx=p0a+p2b
      yy=two*d2b
!cc      a02=(p0a+p2b)*(p0a+p2b)
       a02=(xx)*(xx)
!cc      term1=rsq +  four*d2b*d2b + a02
       term1=rsq +  yy*yy + a02
      term2=rsq + a02
      qxx=1.0D0 / dsqrt(term1)    -      1.0D0 / dsqrt(term2)
      qxx=qxx * half
      twoelcl(2)=qq + qxx
!*******************************

       x1=r+yy
       x2=r-yy
! ( s s | o o)      
!cc      term1=(r + two * d2b)*(r + two * d2b)  +  a02
      term1=(x1)*(x1)  +  a02
      term2=rsq + a02
      term3=(x2)*(x2) + a02
      qoo = pt25 / dsqrt(term1)  -  half/dsqrt(term2)  + pt25 / dsqrt(term3)
      twoelcl(3)=qq+qoo

!  (s s | s o)
       x1=p0a+p1b
       x2=r+d1b
       x3=r-d1b
      a01=(x1)*(x1)
      term1=(x2)*(x2) + a01
      term2=(x3)*(x3) + a01
      quo=one/dsqrt(term1)   -   one/dsqrt(term2)
      quo=half*quo
      twoelcl(4)=quo

! end of hydrogen-heavy atom section+++++++++++++++++++++++



!   heavy atom - hydrogen

      elseif(ni==4.and.nj==1)then
!
!     there are four such integrals as follows (opposite of ones from above)
!     (s s| s s) = [q,q]
!     (p p | s s) = [q,q] + [Qpp , q]
!     (o o | s s) = [q,q] + [oo,q]
!     (s o | s s) = [uo,q]


!     first the (s s | s s) integra
       x1=p0a+p0b
      a00=(x1)*(x1)
      qq=one/dsqrt(rsq+a00)
      twoelcl(1)=qq
!******************************

!   now the (p p | s s ) integral
      x1=p0b+p2a
      x2=two * d2a
      a02=(x1)*(x1)
      term1=rsq +  x2*x2 + a02
      term2=rsq + a02
      qxx=1.0D0 / dsqrt(term1)    -      1.0D0 / dsqrt(term2)
      qxx=qxx * half
      twoelcl(5)=qq + qxx
!******************************


! ( o o | s s)   
       x3=r+x2   
       x4=r-x2
      term1=(x3)*(x3)  +  a02
      term2=rsq + a02
      term3=(x4)*(x4) + a02
      qoo = pt25 / dsqrt(term1)  -  half/dsqrt(term2) &
     + pt25 / dsqrt(term3)
      twoelcl(6)=qq+qoo

!  (s o | s s)
       x1=p0b+p1a
       x2=r + d1a
       x3=r-d1a
      a01=(x1)*(x1)
      term1=(x2)*(x2) + a01
      term2=(x3)*(x3)  + a01
      quo=one/dsqrt(term1)   -   one/dsqrt(term2)
      quo= - half*quo
      twoelcl(7)=quo

! end of heavy atom - hydrogen section+++++++++++++++++++++++


!    heavy atom - heavy atom section

      elseif(ni==4.and.nj==4)then

!
!     there are 22 :(  integrals as follows
!     (s s| s s) = [q,q]
!     (s s | p p) = [q,q] + [q,Qpp]
!     (s s | o o) = [q,q] + [q,oo]
!     (s s | o s) = [q,uo]
      a00=(p0a+p0b)*(p0a+p0b)
      a01=(p0a+p1b)*(p0a+p1b)
      a10=(p1a+p0b)*(p1a+p0b)
      a11=(p1a+p1b)*(p1a+p1b)
      a12=(p1a+p2b)*(p1a+p2b)
      a21=(p2a+p1b)*(p2a+p1b)
      a02=(p0a+p2b)*(p0a+p2b)
      a20=(p2a+p0b)*(p2a+p0b)
      a22=(p2a+p2b)*(p2a+p2b)

!     (s s | s s)
      qq=one/dsqrt(rsq+a00)
      twoelcl(1)=qq

!  (s s | p p)
       x1=two*d2b

      term1=rsq +  x1*x1 + a02
      term2=rsq + a02
      qxx=1.0D0 / dsqrt(term1)    -      1.0D0 / dsqrt(term2)
      qxx=qxx * half
      twoelcl(2)=qq + qxx

! ( s s | o o)
      x2=r+x1
      x3=r-x1
      term1=(x2)*(x2)  +  a02
      term2=rsq + a02
      term3=(x3)*(x3) + a02
      qoo = pt25 / dsqrt(term1)  -  half/dsqrt(term2) &
     + pt25 / dsqrt(term3)
      twoelcl(3)=qq+qoo

! ( s s | s o)

      x1=r+d1b
      x2=r-d1b
      term1=(x1)*(x1) + a01
      term2=(x2)*(x2)  + a01
      quo=one/dsqrt(term1)   -   one/dsqrt(term2)
      quo=half*quo
      twoelcl(4)=quo

! ( p p | s s)

       x1=two*d2a
      term1=rsq +  x1*x1 + a20
      term2=rsq + a20
      qxx=1.0D0 / dsqrt(term1)    -      1.0D0 / dsqrt(term2)
      qxx=qxx * half
      twoelcl(5)=qq + qxx

! ( o o | s s)
       x2=r+x1
       x3=r-x1
      term1=(x2)*(x2)  +  a20
      term2=rsq + a20
      term3=(x3)*(x3) + a20
      qoo = pt25 / dsqrt(term1)  -  half/dsqrt(term2) &
      + pt25 / dsqrt(term3)
      twoelcl(6)=qq+qoo

! ( s o | s s)

       x4=r+d1a     
       x5=r-d1a
      term1=(x4)*(x4) + a10
      term2=(x5)*(x5)  + a10
      quo=one/dsqrt(term1)   -   one/dsqrt(term2)
      quo=-half*quo
      twoelcl(7)=quo


! (s o | p p)
      x1=two*d2b
      term1=(x4)*(x4) + x1*x1 + a12
      term2=(x5)*(x5) + x1*x1 + a12
      term3=(x4)*(x4) + a12
      term4=(x5)*(x5) + a12
      uopp= - pt25 / dsqrt(term1)  +  pt25 / dsqrt(term2) &
     + pt25 / dsqrt(term3) -  pt25 / dsqrt(term4)
      twoelcl(8)=twoelcl(7) + uopp

! (p p | s o)
     x1=r + d1b
     x2=r-d1b
     x3=two *d2a
      term1=(x1)*(x1) + x3*x3 + a21
      term2=(x2)*(x2) + x3*x3 + a21
      term3=(x1)*(x1) + a21
      term4=(x2)*(x2) + a21
      uopp= - pt25 / dsqrt(term1)  +  pt25 / dsqrt(term2) &
     + pt25 / dsqrt(term3) -  pt25 / dsqrt(term4)
      twoelcl(9)=twoelcl(4) - uopp

! (s o | o o)
      x1=r + d1a
      x2=r - d1a
      x3=two*d2b
      yy=x1-x3
      term1=(yy)*(yy) + a12
      yy=x2-x3
      term2=(yy)*(yy) + a12
      yy=x1 + x3
      term3=(yy)*(yy) + a12
      yy=x2+x3
      term4=(yy)*(yy) + a12
      term5=(x1)*(x1) + a12
      term6=(x2)*(x2) + a12
      uooo=- pt125 / dsqrt(term1) + pt125 / dsqrt(term2) - pt125 / dsqrt(term3) + pt125 / dsqrt(term4) & 
     + pt25 / dsqrt(term5) - pt25 / dsqrt(term6)
      twoelcl(10)= twoelcl(7)  +   uooo

! ( o o | o s)
      x1=r + d1b
      x2=r - d1b
      x3=two*d2a
      yy=x1-x3
      term1=(yy)*(yy) + a21
      yy=x2-x3
      term2=(yy)*(yy) + a21
      yy=x1+x3
      term3=(yy)*(yy) + a21
      yy=x2 +x3
      term4=(yy)*(yy) + a21
      term5=(x1)*(x1) + a21
      term6=(x2)*(x2) + a21
      uooo=- pt125 / dsqrt(term1) + pt125 / dsqrt(term2) &
     - pt125 / dsqrt(term3) + pt125 / dsqrt(term4) &
     + pt25 / dsqrt(term5) - pt25 / dsqrt(term6)
      twoelcl(11)= twoelcl(4)  -   uooo

!    ( s p | s p )
      x1=d1a-d1b
      x2=d1a+d1b
      term1=rsq + (x1)*(x1) + a11
      term2=rsq + (x2)*(x2) + a11
      upup=half / dsqrt(term1) - half / dsqrt(term2)
      twoelcl(12)=upup

!  ( s o | s o )
      x1=r+d1a
      x2=r-d1a
      yy=x1-d1b
      term1=(yy)*(yy) + a11
      yy=x1+d1b
      term2=(yy)*(yy) + a11
      yy=x2-d1b
      term3=(yy)*(yy) + a11
      yy=x2+d1b
      term4=(yy)*(yy) + a11
      uouo=pt25 / dsqrt(term1) - pt25 / dsqrt(term2) &
     -pt25 / dsqrt(term3) + pt25 / dsqrt(term4)
      twoelcl(13)=uouo

! (s p | p o)
      x1=r-d2b
      x2=r+d2b
      x3=d1a-d2b
      x4=d1a+d2b
      term1=(x1)*(x1)+(x3)*(x3) + a12
      term2=(x1)*(x1)+(x4)*(x4) + a12
      term3=(x2)*(x2)+(x3)*(x3) + a12
      term4=(x2)*(x2)+(x4)*(x4) + a12
      uppo=-pt25 / dsqrt(term1) + pt25 / dsqrt(term2) &
     +pt25 / dsqrt(term3) - pt25 / dsqrt(term4)
      twoelcl(14)=uppo

! (o p | s p)
      x1=r-d2a
      x2=r+d2a
      x3=d1b-d2a
      x4=d1b+d2a
      term1=(x1)*(x1)+(x3)*(x3)+a21
      term2=(x1)*(x1)+(x4)*(x4)+a21
      term3=(x2)*(x2)+(x3)*(x3)+a21
      term4=(x2)*(x2)+(x4)*(x4)+a21
      uppo=-pt25 / dsqrt(term1) + pt25 / dsqrt(term2) &
     +pt25 / dsqrt(term3) - pt25 / dsqrt(term4)
      twoelcl(15)=-uppo

! ( p o | p o)
      term1=(r+d2a-d2b)*(r+d2a-d2b) + (d2a - d2b)*(d2a - d2b) + a22
      term2=(r+d2a-d2b)*(r+d2a-d2b) + (d2a + d2b)*(d2a + d2b) + a22
      term3=(r+d2a+d2b)*(r+d2a+d2b) + (d2a - d2b)*(d2a - d2b) + a22
      term4=(r+d2a+d2b)*(r+d2a+d2b) + (d2a + d2b)*(d2a + d2b) + a22
      term5=(r-d2a-d2b)*(r-d2a-d2b) + (d2a - d2b)*(d2a - d2b) + a22
      term6=(r-d2a-d2b)*(r-d2a-d2b) + (d2a + d2b)*(d2a + d2b) + a22
      term7=(r-d2a+d2b)*(r-d2a+d2b) + (d2a - d2b)*(d2a - d2b) + a22
      term8=(r-d2a+d2b)*(r-d2a+d2b) + (d2a + d2b)*(d2a + d2b) + a22
      popo= one / dsqrt(term1) - one / dsqrt(term2) &
           -one / dsqrt(term3) + one / dsqrt(term4) &
           -one / dsqrt(term5) + one / dsqrt(term6) &
           +one / dsqrt(term7) - one / dsqrt(term8)
      twoelcl(16)=pt125*popo

! ( p p | p p)
      term1=rsq + four * (d2a - d2b)*(d2a - d2b) + a22
      term2=rsq + four * (d2a + d2b)*(d2a + d2b) + a22
      term3=rsq + four * (d2a)*d2a+ a22
      term4=rsq + four * (d2b)*d2b + a22
      term5=rsq + a22
      pppp=pt125 / dsqrt(term1) + pt125 / dsqrt(term2) &
      -pt25 / dsqrt(term3) - pt25 / dsqrt(term4) &
      +pt25 / dsqrt(term5)
      twoelcl(17)=twoelcl(2) + twoelcl(5) - twoelcl(1) &
      + pppp

! ( p p | p* p*)
      term1=rsq + four *  (d2a*d2a + d2b*d2b) + a22
      term2=rsq + four * d2a*d2a + a22
      term3=rsq + four * d2b*d2b + a22
      term4=rsq + a22
      ppp2p2=pt25 / dsqrt(term1) - pt25 / dsqrt(term2) &
     -pt25 / dsqrt(term3) + pt25 / dsqrt(term4)
      twoelcl(18)=twoelcl(2) + twoelcl(5) - twoelcl(1) &
     +ppp2p2

! ( p p | o o)
      term1=(r - two * d2b)*(r - two * d2b) + four * d2a*d2a + a22
      term2=(r + two * d2b)*(r + two * d2b) + four * d2a*d2a + a22
      term3=(r - two * d2b)*(r - two * d2b) + a22
      term4=(r + two * d2b)*(r + two * d2b) + a22
      term5=rsq + four * d2a*d2a + a22
      term6=rsq + a22
      ppoo=pt125/dsqrt(term1) + pt125 / dsqrt(term2) &
     - pt125/dsqrt(term3) - pt125 / dsqrt(term4) &
     - pt25 / dsqrt(term5) + pt25 / dsqrt(term6)
      twoelcl(19)=twoelcl(3) + twoelcl(5) - twoelcl(1) + ppoo

! ( o o | p p)
      term1=(r - two * d2a)*(r - two * d2a) + four * d2b*d2b + a22
      term2=(r + two * d2a)*(r + two * d2a) + four * d2b*d2b + a22
      term3=(r - two * d2a)*(r - two * d2a) + a22
      term4=(r + two * d2a)*(r + two * d2a) + a22
      term5=rsq + four * d2b*d2b + a22
      term6=rsq + a22
      ppoo=pt125/dsqrt(term1) + pt125 / dsqrt(term2) &
     - pt125/dsqrt(term3) - pt125 / dsqrt(term4) &
     - pt25 / dsqrt(term5) + pt25 / dsqrt(term6)
      twoelcl(20)=twoelcl(6) + twoelcl(2) - twoelcl(1) + ppoo

!  ( o o | o o)
      term1=(r + two * (d2a - d2b))*(r + two * (d2a - d2b)) + a22
      term2=(r + two * (d2a + d2b))*(r + two * (d2a + d2b)) + a22
      term3=(r - two * (d2a + d2b))*(r - two * (d2a + d2b)) + a22
      term4=(r - two * (d2a - d2b))*(r - two * (d2a - d2b)) + a22
      term5=(r + two * d2a)*(r + two * d2a) + a22
      term6=(r - two * d2a)*(r - two * d2a) + a22
      term7=(r + two * d2b)*(r + two * d2b) + a22
      term8=(r - two * d2b)*(r - two * d2b) + a22
      term9=rsq + a22
      oooo=pt06 / dsqrt(term1) + pt06 / dsqrt(term2) &
     +pt06 / dsqrt(term3) + pt06 / dsqrt(term4) &
     -pt125 / dsqrt(term5) - pt125 / dsqrt(term6) &
     -pt125 / dsqrt(term7) - pt125 / dsqrt(term8) &
     +pt25 / dsqrt(term9)
      twoelcl(21)=twoelcl(3) + twoelcl(6) - twoelcl(1) + oooo

!     (p p* | p p*)
!     to maintain rotational invariance, this 
!     is computed as half * [ (pp|pp) - (pp|p*p*)
!
      twoelcl(22)=half * (twoelcl(17) - twoelcl(18))
     
      

      
         


      

      end if
      
!            do 10 ii=1,22
!      twoelcl(ii)=autoev*twoelcl(ii)
!      write(*,*)ii,twoelcl(ii)
! 10   continue

       twoelcl=twoelcl*autoev


      end subroutine local
      
         
      
      
      
         
      
















      

      



