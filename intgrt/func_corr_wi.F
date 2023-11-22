         subroutine func_corr_wi(func)

c The Wilson-Ivanov Functional - Only for closed shell systems

      implicit none
#include <aces.h>
#include <dftfunc.com>

      M_REAL
     &    func

c Parameters in the functional.
c a,b,c,d,k - coefficients
c u - power in numerator
c w - power in denamonator

      M_REAL
     & a,b,c,d,k,u,w,ff,fd

      M_REAL
     & cons83,cons53,term1,term2

      parameter(cons83=8.d0/3.d0,cons53=5.d0/3.d0,
     &    a=-0.06520d0, b=0.00070d0, k=-0.0010d0,
     &    d=0.0010d0,
     &    c=0.210d0,
     &    u=7.d0/2.d0,
     &    w=5.d0) 
         
        if (ro.ge.thresh) then
          ff=gro2/(ro**cons83)
          fd=1.d0+c*rs+(d*(gro**u)/(ro**w))
          term1=a*ro/fd
          term2=b*(ff*ro)*dexp(k*ff)/fd
          func=term1+term2
        else
          func=0.d0
        end if

        return
        end
