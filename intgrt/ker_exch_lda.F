c LDA exchange kernel
c
c                              1/3    -2/3
c      fxa = - (4/3) (1/3) Cx 2   (roa   )
c
c                           1/3
c      Cx=(9/8) alpha (3/pi)

         subroutine ker_exch_lda(kerna,kernb)

      implicit none
#include <aces.h>
#include <iuhf.com>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    kerna,kernb

      M_REAL
     &    con13,con43,con98,alpha,conlda,pi
      parameter (con13=1.d0/3.d0,con43=4.d0/3.d0,con98=9.d0/8.d0)
      parameter (alpha=2.d0/3.d0)

      callstack_curr='KER_EXCH_LDA'

      pi = acos(-1.d0)
      conlda=alpha*(2.d0**con13)*con98*((3.d0/pi)**con13)

      if (roa.ge.thresh) then
         kerna=-con43*con13*conlda*(roa**(-alpha))
      else
         kerna=0.d0
      end if

      if (iuhf.eq.0) then
         kernb=kerna
      else
         if (rob.ge.thresh) then 
           kernb=-con43*con13*conlda*(rob**(-alpha))
         else
           kernb=0.d0
         end if
      end if

      return
      end
