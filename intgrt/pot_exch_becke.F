c     Becke exchange potential
c
c     by Stan Ivanov, January 1999

      subroutine pot_exch_becke(pota,grcompa,potb,grcompb)

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    pota,grcompa,potb,grcompb

      M_REAL
     &    potbeckea,potbeckeb,
     &    groa, grob

      M_REAL
     &    x,x2,f,g,asinh,
     &    con43,con13,beta

        parameter (con43=4.d0/3.d0,
     &            con13=1.d0/3.d0)

        parameter(beta=0.0042d0)

      callstack_curr='POT_EXCH_BECKE'

      call pot_exch_lda(pota,potb)

c Becke exchange functional has two pieces.
c We find the second piece that involves the
c the gradient corrections in beckepot.
c The first piece is the LDA exchange.
c The second one is the Becke correction f.
c potbecke is df/drho piece and
c grcomp is df/d|gradrho| piece.

      grcompa=0.d0
      grcompb=0.d0

       if (roa.ge.thresh) then
         groa=dsqrt(gro2a)
         x=groa/(roa**con43)
         x2=x*x

         asinh=dlog(x+dsqrt(x2+1.d0))
         g=(1.d0+6.d0*beta*x*asinh)
         f=dsqrt(1.d0+(x2))

         potbeckea=-beta*
     &     ((con43*(roa**con13)*x2*(beta*6.d0*x2-f))/
     &      (g*g*f))

         grcompa=-beta*((-2.d0*x*(-f-
     &     3.d0*beta*x*f*asinh+3.d0*beta*x2))/(f*g*g))

         pota=pota+potbeckea
       end if


      if (iuhf.eq.0) then
        potb=pota
        grcompb=grcompa
      else
        if (rob.ge.thresh) then
          grob=dsqrt(gro2b)
          x=grob/(rob**con43)
          x2=x*x

          asinh=dlog(x+dsqrt(x2+1.d0))
          g=(1.d0+6.d0*beta*x*asinh)
          f=dsqrt(1.d0+(x2))

          potbeckeb=-beta*
     &      ((con43*(rob**con13)*x2*(beta*6.d0*x2-f))/
     &      (g*g*f))

          grcompb=-beta*((-2.d0*x*(-f-
     &      3.d0*beta*x*f*asinh+3.d0*beta*x2))/(f*g*g))

          potb=potb+potbeckeb
        end if

      endif
      return
      end
