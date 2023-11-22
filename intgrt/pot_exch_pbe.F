c       PBE exchange potential  
c                                          
c       by Stan Ivanov, January 1999              
c   
c       Ref.: Phys. Rev. Lett. 77, 3865 (1996), ibid 78, 1396 (1997). 

         subroutine pot_exch_pbe(pota,grcompa,potb,grcompb)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <iuhf.com>
#include <callstack.com>

      M_REAL
     &    pota,grcompa,potb,grcompb

      M_REAL
     &    groa, grob


      M_REAL
     &    con13,con43,con83,pi

      parameter (con13=1.d0/3.d0,con43=4.d0/3.d0,con83=8.d0/3.d0)

      M_REAL
     &    s2,fx,dfxds,mu,k,conlda,cons

      parameter(mu=0.2195149727645171d0,k=0.8040d0)

      fx(s2)=1.d0+k-k/(1.d0+mu*s2/k)
      dfxds(s2)=2.d0*mu*dsqrt(s2)/((1.d0+mu*s2/k)*(1.d0+mu*s2/k))

      callstack_curr='POT_EXCH_PBE'

      pi = acos(-1.d0)
      conlda=-(81.0d0/(32.0d0*pi))**con13
      cons=(48.0d0*pi*pi)**(-con13)


      if (roa.ge.thresh) then
         groa=dsqrt(gro2a)
         s2=cons*cons*gro2a/(roa**con83)

         pota=con43*conlda*(roa**con13)*fx(s2)
     &        -conlda*dfxds(s2)*cons*groa*con43*(1.d0/roa)
         grcompa=conlda*dfxds(s2)*cons        
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
             s2=cons*cons*gro2b/(rob**con83)

             potb=con43*conlda*(rob**con13)*fx(s2)
     &        -conlda*dfxds(s2)*cons*grob*con43*(1.d0/rob)

             grcompb=conlda*dfxds(s2)*cons
           else
             potb=0.d0
             grcompb=0.d0
          endif
      end if

      return
      end
