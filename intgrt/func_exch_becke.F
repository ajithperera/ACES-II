c       Becke exchange functional  
c                                          
c       by Stan Ivanov, January 1999              

      subroutine func_exch_becke(func)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <iuhf.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &    con13,con32,con43,con83,beta
      parameter (con13=1.d0/3.d0,con32=3.d0/2.d0,con43=4.d0/3.d0,
     &    con83=8.d0/3.d0)
      parameter (beta=0.0042d0)

      M_REAL
     &    asinh,g,x,x2,xa2,xb2,xa,xb,
     &    fbeckea,fbeckeb

      asinh(x)=dlog(x+dsqrt(x*x+1.d0))
      g(x,x2)=beta*x2/(1.d0+6.d0*beta*x*asinh(x))

      callstack_curr='FUNC_EXCH_BECKE'

      call func_exch_lda(func)

      if (roa.ge.thresh) then
        xa2=gro2a*roa**(-con83)
        xa=dsqrt(xa2)
        fbeckea=-((roa**con43)*g(xa,xa2))
      else
        fbeckea=0.d0
      endif

      if (iuhf.eq.0) then
        fbeckeb=fbeckea
      else 
          if (rob.ge.thresh) then
             xb2=gro2b*rob**(-con83)
             xb=dsqrt(xb2)
             fbeckeb=-((rob**con43)*g(xb,xb2))
          else
             fbeckeb=0.d0 
          endif
      end if

      func=func+fbeckea+fbeckeb

      return
      end
