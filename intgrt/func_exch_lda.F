c LDA exchange energy.
c
c                1/3    4/3     4/3
c      E = - Cx 2   (roa   + rob   )
c   
c LDA exchange potential.
c
c                       1/3    1/3
c      Va = - (4/3) Cx 2   (roa   )
c
c                           1/3
c      Cx=(9/8) alpha (3/pi)
c
c R. G. Parr and W. Yang, Density-Functional Theory of Atoms and Molecules,
c             Oxford Science Publications, 1989.  p. 156
c
c Modified by Stan Ivanov

      subroutine func_exch_lda(func)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <iuhf.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &    con13,con43,con98,alpha,conlda,funa,funb,pi
      parameter (con13=1.d0/3.d0,con43=4.d0/3.d0,con98=9.d0/8.d0)
      parameter (alpha=2.d0/3.d0)

      callstack_curr='FUNC_EXCH_LDA'

      pi = acos(-1.d0)
      conlda=alpha*(2.d0**con13)*con98*((3.d0/pi)**con13)

      funa=0.d0
      if (roa.ge.thresh) funa=-conlda*roa**con43

      funb=funa
      if (iuhf.ne.0) then
         funb=0.d0
         if (rob.ge.thresh) funb=-conlda*rob**con43
      end if

      func=funa+funb

      return
      end

