      function dldaf(roa,rob,grdaa,grdab,grdbb)
      
c Determine the value of the LDA functional.

      implicit none

c *** FUNCTION ARGUMENTS *********************************************

      double precision
     &    dldaf,roa,rob,grdaa,grdab,grdbb

c *** LOCAL VARIABLES ************************************************

      double precision
     &    ro,x,zeta,g,ea,ep,ef,
     &    xx,t1,t2,t3,t4,t5,e,a,b,c,q,x0,pi

      double precision
     &    aa,ba,ca,x0a,qa,ap,bp,cp,x0p,qp,af,bf,cf,x0f,qf

c ********************************************************************
c ********************************************************************

      xx(x   ,b,c)=x*x+b*x+c
      t1(x   ,b,c)=log(x*x/xx(x,b,c))
      t2(x   ,b,q)=2*b/q*atan(q/(2*x+b))
      t3(  x0,b,c)=b*x0/xx(x0,b,c)
      t4(x,x0,b,c)=log((x-x0)**2/xx(x,b,c))
      t5(x,x0,b,q)=2*(2*x0+b)/q*atan(q/(2*x+b))
      e(a,b,c,x,x0,q)=a*(t1(x,b,c)+t2(x,b,q)-t3(x0,b,c)*
     &    (t4(x,x0,b,c)+t5(x,x0,b,q)))


      pi=4.0*atan(1.0)
      aa=-1.0/(3.0*pi*pi)
      ba=1.13107
      ca=13.0045
      x0a=-0.00475840
      qa=sqrt(4.0*ca-ba*ba)

      ap=0.0621814
      bp=3.72744
      cp=12.9352
      x0p=-0.10498
      qp=sqrt(4.0*cp-bp*bp)

      af=0.5*ap
      bf=7.06042
      cf=18.0578
      x0f=-0.32500
      qf=sqrt(4.0*cf-bf*bf)

      ro=roa+rob
      if (ro.lt.1.0d-20) then
        dldaf=0.0
        return
      endif

      x=(3.0/(4.0*pi*ro))**(1.0/6.0)
      zeta=(roa-rob)/ro
      g=9.0/8.0*((1.0+zeta)**(4.0/3.0)+(1.0-zeta)**(4.0/3.0)-
     &    2.0)

      ea=e(aa,ba,ca,x,x0a,qa)
      ep=e(ap,bp,cp,x,x0p,qp)
      ef=e(af,bf,cf,x,x0f,qf)

C      dldaf=1.0+(4.0/(9.0*(2.0**(1.0/3.0)-1.0))*
C     &    (ef-ep)/ea-1.0)*zeta**4.0
C Ajith 05/29/96 DEC Alpha do not like this
C
      dldaf=1.0+(4.0/(9.0*(2.0**(1.0/3.0)-1.0))*
     &      (ef-ep)/ea-1.0)*zeta**4
C
c The equations are in rydberg
      dldaf=0.5*ro*(ep+ea*g*dldaf)
      
      return
      end
