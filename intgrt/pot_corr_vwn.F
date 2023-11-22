
      subroutine pot_corr_vwn(pota,potb)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    pota,potb

      M_REAL
     &    con43,con13,con98,con49,con32,pi,
     &    x0p,x0f,x0a,cp,cf,ca,bp,bf,ba,ap,af,aa

      parameter(con43=4.d0/3.d0,con13=1.d0/3.d0,con98=9.d0/8.d0,
     &    con49=4.d0/9.d0,con32=3.d0/2.d0,
     &    x0p=-0.10498d0, ap=0.03109070d0,   bp=3.72744d0, cp=12.9352d0,
     &    x0f=-0.325d0,   af=0.5d0*ap,       bf=7.06042d0, cf=18.0578d0,
     &    x0a=-0.0047584d0,ba=1.13107d0,
     &    ca=13.0045d0)

      M_REAL
     &    qi,b,c,xxi,x,ei,a,x0,q,xx,xx0,delei,qp,qf,qa,xxp,xxf,xxa,
     &    xx0p,xx0f,xx0a,ep,ef,ea,delep,delef,delea,con,
     &    g,h,delg,delh,delx,delzeta,delzetb,deleca,delecb,fact1,fact2,
     &    zeta4,ec

      qi(b,c)=dsqrt(4.d0*c-b*b)
      xxi(x,b,c)=x*x+b*x+c
      ei(a,x,x0,b,q,xx,xx0)=a*(dlog(x*x/xx) +
     &    (b+b)/q*datan(q/(x+x+b)) -
     &    b*x0/xx0* (dlog((x-x0)**2/xx)+
     &    (2.d0*(b+x0+x0)/q)*datan(q/(x+x+b))  )  )
      delei(a,x,x0,b,q,xx,xx0)=a*(2/x-(2*x+b)/xx-4*b/((2*x+b)**2+q**2)-
     &    b*x0/xx0*(2.d0/(x-x0)-(2*x+b)/xx-
     &    4.d0*(2.d0*x0+b)/((2*x+b)**2+q**2)))

      callstack_curr='POT_CORR_VWN'

      pi = acos(-1.d0)
      aa = -1.d0/(6.d0*pi*pi)

      if (ro.ge.thresh) then

        qp=qi(bp,cp)
        qf=qi(bf,cf)
        qa=qi(ba,ca)
        x=sqrt(rs)
        xxp=xxi(x,bp,cp)
        xxf=xxi(x,bf,cf)
        xxa=xxi(x,ba,ca)
        xx0p=xxi(x0p,bp,cp)
        xx0f=xxi(x0f,bf,cf)
        xx0a=xxi(x0a,ba,ca)
        ep=ei(ap,x,x0p,bp,qp,xxp,xx0p)
        ef=ei(af,x,x0f,bf,qf,xxf,xx0f)
        ea=ei(aa,x,x0a,ba,qa,xxa,xx0a)
        delep=delei(ap,x,x0p,bp,qp,xxp,xx0p)
        delef=delei(af,x,x0f,bf,qf,xxf,xx0f)
        delea=delei(aa,x,x0a,ba,qa,xxa,xx0a)

        con=con49/(2.d0**con13-1.d0)
        g=con98*((1.d0+zeta)**con43+(1.d0-zeta)**con43-2.d0)
        delg=con32*((1.d0+zeta)**con13-(1.d0-zeta)**con13)
        h=con*(ef-ep)/ea-1.d0
        delh=con/ea*(delef-delep-(ef-ep)*delea/ea)
        delx=-(x/ro)/6.d0
        delzeta=(1.d0-zeta)/ro
        delzetb=-(1.d0+zeta)/ro
        zeta4=zeta**4

        fact1=delx*(delep+delea*g*(1.d0+h*zeta4)+ea*g*delh*zeta4)
        fact2=ea*(delg*(1.d0+h*zeta4)+4.d0*g*h*zeta**3)
        deleca=fact1+delzeta*fact2
        delecb=fact1+delzetb*fact2
        ec=ep+ea*g*(1.d0+h*zeta4)

        pota=ec+ro*deleca
        potb=ec+ro*delecb
      else
        pota=0.d0
        potb=0.d0
      end if

      return
      end
