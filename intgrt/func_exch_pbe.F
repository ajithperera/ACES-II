c       PBE exchange functional  
c                                          
c       by Stan Ivanov, January 1999              
c   
c       Ref.: Phys. Rev. Lett. 77, 3865 (1996), ibid 78, 1396 (1997). 

      subroutine func_exch_pbe(func)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <iuhf.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &    con13,con43,con83

      parameter (con13=1.d0/3.d0,con43=4.d0/3.d0,con83=8.d0/3.d0)

      M_REAL
     &    fpbea,fpbeb,s2,fx,mu,k,conlda,cons,pi

      parameter(mu=0.2195149727645171d0,k=0.8040d0)

      fx(s2)=1.d0+k-k/(1.d0+mu*s2/k)

      callstack_curr='FUNC_EXCH_PBE'

      pi = acos(-1.d0)
      conlda=-(81.0d0/(32.0d0*pi))**con13
      cons=(48.0d0*pi*pi)**(-con13)

      if (roa.ge.thresh) then
         s2=cons*cons*gro2a/(roa**con83)
         fpbea=conlda*(roa**con43)*fx(s2)
      else
         fpbea=0.d0 
      endif

      if (iuhf.eq.0) then
          fpbeb=fpbea
      else 
          if (rob.ge.thresh) then
             s2=cons*cons*gro2b/(rob**con83)
             fpbeb=conlda*(rob**con43)*fx(s2)
          else
             fpbeb=0.d0
          endif
      end if

      func=fpbea+fpbeb

      return
      end
