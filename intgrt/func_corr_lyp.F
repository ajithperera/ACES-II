      subroutine func_corr_lyp(func)

c    Lee-Yang-Parr correlation potential
c    Written by S. Ivanov, February 1999
c
c    Chem. Phys. Lett. 157, 200 (1989)

 
      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &   func

        M_REAL
     &    con13,con23,con83,con113,pi,
     &    ac,bc,cc,dc,cf

      parameter (con13=1.d0/3.d0,con23=2.d0/3.d0,
     &    con83=8.d0/3.d0, con113=11.0d0/3.d0)

      parameter (ac=0.04918d0,bc=0.132d0,cc=0.2533d0,dc=0.349d0)

      M_REAL
     &    ro13m,denom,t1den,t1,omega,delta,fcomm,f,hcomm,h,
     &    roa_rob

      callstack_curr='FUNC_CORR_LYP'

      pi = acos(-1.d0)

      if (ro.ge.thresh) then

        cf=0.30d0*(3.0d0*pi*pi)**con23

        ro13m=ro**(-con13)
        denom=1.d0+dc*ro13m
        t1den=ro*denom
        t1=-ac*4.d0/t1den

        omega=dexp(-cc*ro13m)/(denom*ro**con113)
        delta=ro13m*(cc+dc/denom)

        fcomm=(roa*gro2a+rob*gro2b)*roinv/9.0d0
        f=(2.d0**con113)*cf*(roa**con83+rob**con83)+
     &    (47.0d0-7.0d0*delta)*gro2/18.0d0-
     &    (2.50d0-delta/18.0d0)*(gro2a+gro2b)-
     &    (delta-11.0d0)*fcomm

        hcomm=con23*ro*ro
        h=-hcomm*gro2+
     &         (hcomm-roa*roa)*gro2b+(hcomm-rob*rob)*gro2a

        roa_rob=roa*rob
        func=roa_rob*t1-ac*bc*omega*(roa_rob*f+h)

      else
        func=0.d0 
      end if

      return
      end
