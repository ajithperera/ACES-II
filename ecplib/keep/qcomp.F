      function qcomp(n,l)
c     ----- computes q(n,l)                                -----
c     ----- scaled by exp(-t) to prevent overflows         -----
c     ----- arguments are alpha, xk, and t=xk**2/(4*alpha) -----
c     ----- no restriction on the magnitude of t           -----
c     ----- increase dfac array to raise n, l restrictions -----
      implicit double precision (a-h,o-z)

      common /pifac/ pi, sqpi2
      common /dfac/ dfac(29)
      common/qstore/dum1(81),alpha,xk,t
      dimension tmin(9)

      data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,
     1   18.0d0,15.0d0/
      data am1,   a0,   accpow, accasy, a1,   a2,   a4
     1    /-1.0d0,0.0d0,1.0d-14,1.0d-10,1.0d0,2.0d0,4.0d0/

      if(mod(n+l,2).ne.0.or.n.le.l) go to 30
c     ----- use alternating series (n+l.le.22.and.l.le.10) -----

      call cortab
      if(l.eq.0) xkp=a1
      if(l.ne.0) xkp=(xk/(alpha+alpha))**l
      prefac=sqpi2*xkp*dfac(n+l+1)/
     1   ((alpha+alpha)**((n-l)/2)*dsqrt(alpha+alpha)*dfac(l+l+3))
      num=l-n+2
      xden=(l+l+3)
      term=a1
      sum=term
      xc=am1
 10   if(num.eq.0) go to 20
      fnum=num
      term=term*fnum*t/(xden*xc)
      xc=xc+am1
      sum=sum+term
      num=num+2
      xden=xden+a2
      go to 10
   20 qcomp=prefac*sum
      return
   30 if(t.lt.tmin(min0(n,8)+1)) go to 60
c     ----- use asymptotic series (arbitrary n,l) -----
      xkp=(xk/(alpha+alpha))**(n-2)
      prefac=xkp*sqpi2/((alpha+alpha)*dsqrt(alpha+alpha))
      sum=a1
      term=a1
      fac1=(l-n+2)
      fac2=(1-l-n)
      xc=a1
   40 term=term*fac1*fac2/(a4*xc*t)
      if(term.eq.a0) go to 50
      sum=sum+term
      if(dabs(term/sum).lt.accasy) go to 50
      fac1=fac1+a2
      fac2=fac2+a2
      xc=xc+a1
      go to 40
   50 qcomp=prefac*sum
      return
c     ----- use power series (n+l.le.22.and.l.le.10) -----
   60 if(l.eq.0) xkp=a1
      if(l.ne.0) xkp=(xk/(alpha+alpha))**l
      prefac=dexp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
      if(mod(n+l,2).eq.0) prefac=prefac*sqpi2/dsqrt(alpha+alpha)
      xnum=(l+n-1)
      xden=(l+l+1)
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=a0
   70 xnum=xnum+a2
      xden=xden+a2
      xj=xj+a1
      term=term*t*xnum/(xj*xden)
      sum=sum+term
      if((term/sum).gt.accpow) go to 70
      qcomp=prefac*sum
      return
      end
      subroutine cortab
      implicit double precision (a-h,o-z)
      common /pifac/ pi, sqpi2
      common /dfac/ dfac(29)
c*    common /dfac/ dfac(23)
c*    common/fact/fac(13),fprod(7,7)
      common/fact/fac(17),fprod(9,9)
c*    common /zlmtab/ zlm(130)
      common /zlmtab/ zlm(377)
      data a0, a1s2, a1, a2, a3 /0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0/
c     ----- this value of pi may have to be changed if the     -----
c     ----- compiler will not accept extra significant figures -----
      pi=3.141592653589793238462643383279502884197169399d0
      sqpi2=sqrt(a1s2*pi)
      rt3=sqrt(a3)
      rt6=sqrt(6.0d0)
      rt5s2=sqrt(2.5d0)
      rt3s2=sqrt(1.5d0)
      dfac(1)=a1
      dfac(2)=a1
      fi=a0
c*    do 10 i=1,21
      do 10 i=1,27
      fi=fi+a1
   10 dfac(i+2)=fi*dfac(i)
      fac(1)=a1
      fi=a0
c*    do 20 i=1,12
      do 20 i=1,16
      fi=fi+a1
   20 fac(i+1)=fi*fac(i)
c*    do 30 l1=1,7
      do 30 l1=1,9
      do 30 k1=1,l1
   30 fprod(k1,l1)=fac(l1+(k1-1))/(fac(k1)*fac(l1-(k1-1)))
      Return 
      End

     
