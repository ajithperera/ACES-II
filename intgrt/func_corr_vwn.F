c Vosko, Wilk, Nusair correlation energy.
c
c Written by S. Beck; Modified by S. Ivanov
c
c
c    [4.4] defines e(p), e(f), and e(a)
c    [2.4] Ec(r,zeta)     = e(p) + del ec(r,zeta)
c    [4.7] del ec(r,zeta) = ...
c          del ec(r,1)/del eRPA(r,1) ~ 1   (see table 5)
c    [3.2] del eRPA(r,zeta) = ...
c          alphaRPA == e(a)
c    [2.3] f(zeta) =
c          f''(0) = 8/9 [2**4/3 - 2]**-1
c          betaRPA(r) = beta1(r)
c    [4.8] beta1(r) =
c          del ec(r,1) = e(f) - e(p)   (see table 5)
c
c Note that the energy units in the paper are in Rydberts.  To switch
c to Hartrees, the values of A are all divided by 2.
c
c S. H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys., 58, 1200, 1980.
c
      
           subroutine func_corr_vwn(func)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &    con43,con13,con98,pi,
     &    x0p,x0f,x0a,cp,cf,ca,bp,bf,ba,ap,af,aa

      parameter(con43=4.d0/3.d0,con13=1.d0/3.d0,con98=9.d0/8.d0,
     &    x0p=-0.10498d0, ap=0.03109070d0,   bp=3.72744d0, cp=12.9352d0,
     &    x0f=-0.325d0,   af=0.5d0*ap,       bf=7.06042d0, cf=18.0578d0,
     &    x0a=-0.0047584d0,ba=1.13107d0,
     &    ca=13.0045d0)

      M_REAL
     &    qi,xxi,b,c,x,x0,g1,g2,
     &    ei,a,q,xx,xx0,qp,qf,qa,xs,xxp,xxf,
     &    xxa,xx0p,xx0f,xx0a,ep,ef,ea

      qi(b,c)=dsqrt(4.d0*c-b*b)
      xxi(x,b,c)=x*x+b*x+c
      ei(a,x,x0,b,q,xx,xx0)=a*(dlog(x*x/xx) +
     &    (b+b)/q*datan(q/(x+x+b)) -
     &    b*x0/xx0* (dlog((x-x0)**2/xx)+
     &    (2.d0*(b+x0+x0)/q)*datan(q/(x+x+b))  )  )

      callstack_curr='FUNC_CORR_VWN'

      pi = acos(-1.d0)
      aa = -1.d0/(6.d0*pi*pi)

      func=0.d0
      if (ro.ge.thresh) then 
         qp=qi(bp,cp)
         qf=qi(bf,cf)
         qa=qi(ba,ca)
         xs=dsqrt(rs)
         xxp=xxi(xs,bp,cp)
         xxf=xxi(xs,bf,cf)
         xxa=xxi(xs,ba,ca)
         xx0p=xxi(x0p,bp,cp)
         xx0f=xxi(x0f,bf,cf)
         xx0a=xxi(x0a,ba,ca)
         ep=ei(ap,xs,x0p,bp,qp,xxp,xx0p)
         ef=ei(af,xs,x0f,bf,qf,xxf,xx0f)
         ea=ei(aa,xs,x0a,ba,qa,xxa,xx0a)

         g1=4.d0/(9.d0*(2.d0**con13-1.d0))
         g2=con98*((1.d0+zeta)**con43+(1.d0-zeta)**con43-2.d0)
         func=1.d0+(g1*(ef-ep)/ea-1.d0)*zeta**4

         func=ro*(ep+ea*g2*func)
      end if

      return
      end

