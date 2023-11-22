         subroutine func_corr_wl(func)

c The Wilson-Levy Functional - only the functional.
c WL is best as a "tack-on" functional. 

      implicit none
#include <aces.h>
#include <dftfunc.com>

      M_REAL
     &    func

c Parameters in the functional.
c a,b,c,d - coefficients

      M_REAL
     & a,b,c,d

      M_REAL
     & cons43,term1,term2,pol,xa,xb

      parameter(cons43=4.d0/3.d0,
     &    a=-0.74860d0, b=0.060010d0, 
     &    c=3.600730d0,
     &    d=0.90d0)
         
        if (ro.ge.thresh) then
          pol=dsqrt(1.d0-zeta*zeta)
          term1=ro*(a+b*gro/(ro**cons43))

          if (roa.lt.thresh) then
             xa=0.d0
          else
             xa=dsqrt(gro2a)/(roa**cons43)
          end if
          if (rob.lt.thresh) then
             xb=0.d0
          else
             xb=dsqrt(gro2b)/(rob**cons43)
          end if

          term2=c+rs+d*(xa+xb)
          func=pol*term1/term2
        else
          func=0.d0 

        end if

        return
        end
