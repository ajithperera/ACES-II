c       PW91 exchange functional  
c                                          
c       by Stan Ivanov, January 1999              
c   
c     Phys. Rev. B 54, 16 533 (1996).
c     Phys. Rev. B 45, 13 244 (1992).

      subroutine func_exch_pw91(func)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <iuhf.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &    con13,con43,con83,pi

      parameter (con13=1.d0/3.d0,con43=4.d0/3.d0,con83=8.d0/3.d0)

      M_REAL
     &    fpwa,fpwb,s,s2,fx,fx1,conlda,cons,a,a1,a2,a3,a4,b1,t

      parameter(a=7.79560d0,a1=0.196450d0,a2=0.27430d0,a3=0.150840d0,
     &          a4=100.0d0,b1=0.0040d0)

      fx1(s,s2)=1.d0+a1*s*dlog(a*s+dsqrt(1+a*a*s2))

      callstack_curr='FUNC_EXCH_PW91'

      pi = acos(-1.d0)
c       In the code by K. Burke the value of -0.7385588 was used.
c       We shall use more accurate value for conlda.
      conlda=-(81.0d0/(32.0d0*pi))**con13

c      Old one:
c      conlda=-0.7385588*2.d0**con13
      cons=(48.0d0*pi*pi)**(-con13)

      if (roa.ge.thresh) then
         s2=cons*cons*gro2a/(roa**con83)
         s=dsqrt(s2)
         t=fx1(s,s2)
         fx=(t+a2*s2-a3*s2*dexp(-a4*s2))/(t+b1*s2*s2)
         fpwa=conlda*(roa**con43)*fx
      else
         fpwa=0.d0 
      endif

      if (iuhf.eq.0) then
         fpwb=fpwa
      else 
         if (rob.ge.thresh) then
             s2=cons*cons*gro2b/(rob**con83)
             s=dsqrt(s2)
             t=fx1(s,s2)
             fx=(t+a2*s2-a3*s2*dexp(-a4*s2))/(t+b1*s2*s2) 
             fpwb=conlda*(rob**con43)*fx
          else
             fpwb=0.d0
          endif
      end if

      func=fpwa+fpwb

      return
      end
