
      subroutine pot_exch_lda(pota,potb)
      implicit none
#include <aces.h>
#include <iuhf.com>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    pota,potb

      M_REAL
     &    con13,con43,con98,alpha,conlda,pi
      parameter (con13=1.d0/3.d0,con43=4.d0/3.d0,con98=9.d0/8.d0)
      parameter (alpha=2.d0/3.d0)

      callstack_curr='POT_EXCH_LDA'

      pi = acos(-1.d0)

      conlda=alpha*(2.d0**con13)*con98*((3.d0/pi)**con13)

      if (roa.ge.thresh) then
         pota=-con43*conlda*(roa**con13)
      else
         pota=0.d0
      end if

      if (iuhf.eq.0) then
         potb=pota
      else
         if (rob.ge.thresh) then 
           potb=-con43*conlda*(rob**con13)
         else
           potb=0.d0
         end if
      end if

      return
      end
