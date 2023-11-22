c       PW91 exchange potential  
c                                          
c       by Stan Ivanov, January 1999              
c   
c     Phys. Rev. B 54, 16 533 (1996).
c     Phys. Rev. B 45, 13 244 (1992).

         subroutine pot_exch_pw91(pota,grcompa,potb,grcompb)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <iuhf.com>
#include <callstack.com>

      M_REAL
     &    pota,grcompa,potb,grcompb

      M_REAL
     &    con13,con43,con83,pi

      parameter (con13=1.d0/3.d0,con43=4.d0/3.d0,con83=8.d0/3.d0)

      M_REAL
     &    groa,grob,s,s2,
     &    fxsr,fxlog,t,dtds,fx,dfxds,
     &    num,den,dnumds,ddends

      M_REAL
     &    conlda,cons,a,a1,a2,a3,a4,b1

      parameter(a=7.79560d0,a1=0.196450d0,a2=0.27430d0,a3=0.150840d0,
     &          a4=100.0d0,b1=0.0040d0)


c    Formulas

c      t=1.d0+a1*s*dlog(a*s+dsqrt(1+a*a*s2))   
c      dtds=a1*dlog(a*s+dsqrt(1+a*a*s2))+
c     &  a1*s*(1.d0/(a*s+dsqrt(1+a*a*s2)))*(a+(a*a*s/dsqrt(1+a*a*s2)))

      callstack_curr='POT_EXCH_PW91'

      pi = acos(-1.d0)

c       The following value of conlda is more accurate than what was
c       used in the original code by K. Burke.
      conlda=-(81.0d0/(32.0d0*pi))**con13
c       Old one:
c       conlda=-0.7385588*2.d0**con13
      cons=(48.0d0*pi*pi)**(-con13)

      if (roa.ge.thresh) then
         groa=dsqrt(gro2a)
         s=cons*groa/(roa**con43)

         s2=s*s
         fxsr=dsqrt(1+a*a*s2)
         fxlog=a1*dlog(a*s+fxsr)
         t=1.d0+s*fxlog

         num=t+a2*s2-a3*s2*dexp(-a4*s2)
         den=t+b1*s2*s2
         fx=num/den 

         dtds=fxlog+a1*s*(1.d0/(a*s+fxsr))*(a+(a*a*s/fxsr))

         dnumds=dtds+2.d0*s*a2-2.d0*a3*s*dexp(-a4*s2)+
     &           2.d0*a3*a4*s2*s*dexp(-a4*s2)
         ddends=dtds+4.d0*b1*s2*s
         dfxds=(dnumds/den)-(num*ddends)/(den*den)

         pota=con43*conlda*(roa**con13)*fx
     &      -conlda*dfxds*cons*groa*con43*(1.d0/roa)
         grcompa=conlda*dfxds*cons        
      else
         pota=0.d0
         grcompa=0.d0
      endif

      if (iuhf.eq.0) then
         potb=pota
         grcompb=grcompa
      else 
         if (rob.ge.thresh) then
             grob=dsqrt(gro2b)
             s=cons*grob/(rob**con43)
 
             s2=s*s
             fxsr=dsqrt(1+a*a*s2)
             fxlog=a1*dlog(a*s+fxsr)
             t=1.d0+s*fxlog

             num=t+a2*s2-a3*s2*dexp(-a4*s2)
             den=t+b1*s2*s2
             fx=num/den

             dtds=fxlog+a1*s*(1.d0/(a*s+fxsr))*(a+(a*a*s/fxsr))

             dnumds=dtds+2.d0*s*a2-2.d0*a3*s*dexp(-a4*s2)+
     &           2.d0*a3*a4*s2*s*dexp(-a4*s2)
             ddends=dtds+4.d0*b1*s2*s
             dfxds=(dnumds/den)-(num*ddends)/(den*den) 

             potb=con43*conlda*(rob**con13)*fx
     &          -conlda*dfxds*cons*grob*con43*(1.d0/rob)
             grcompb=conlda*dfxds*cons
         else
             potb=0.d0
             grcompb=0.d0
         endif
      end if

      return
      end
