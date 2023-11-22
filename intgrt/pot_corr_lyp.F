      subroutine pot_corr_lyp(pota,grcompa,gracompa,
     &                        potb,grcompb,grbcompb)

c    Lee-Yang-Parr correlation potential
c    Written by S. Ivanov, February 1999
c    For LYP we need tree different quantities:
c    dF/drho(alpha,beta), dF/d|grad_rho|, 
c    and dF/d|grad_rho(alpha,beta)| returned as
c    pot, grcomp, and gr*comp respectively.
c
c    Chem. Phys. Lett. 157, 200 (1989)

 
      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    pota,potb,grcompa,grcompb,gracompa,grbcompb

        M_REAL
     &    con13,con23,con43,con53,con83,con103,con113,pi,
     &    ac,bc,cc,dc,cf

      parameter (con13=1.d0/3.d0,con23=2.d0/3.d0,con43=4.d0/3.d0,
     &    con53=5.d0/3.d0,con83=8.d0/3.d0,
     &    con103=10.d0/3.d0, con113=11.0d0/3.d0)
      parameter (ac=0.04918d0,bc=0.132d0,cc=0.2533d0,dc=0.349d0)
c      parameter (cf=2.8711234d0)

      M_REAL
     &   denom,t1den,t1,t1ro,
     &   omega,omegaro,delta,deltaro,
     &   fcon,fcomm,f,fro,frocomp1,frocomp2,fgro,fgrocomp1,fgrocomp2,
     &   fq1,fq2,fq3,
     &   hcomm,h,hcommro,hro,hgro,
     &   ab_omega,roa_rob,comm,groa,grob,
     &   ro13,rom13,ro43,rom43

      callstack_curr='POT_CORR_LYP'

      pi = acos(-1.d0)

      if (ro.ge.thresh) then

        cf=0.30d0*(3.d0*pi*pi)**con23
        roa_rob=roa*rob
        ro13=ro**con13
        rom13=1.d0/ro13
        ro43=ro**con43
        rom43=1.d0/ro43
        groa=dsqrt(gro2a)
        grob=dsqrt(gro2b)

        denom=1.d0+dc*rom13

        t1den=ro*denom
        t1=-ac*4.d0/t1den
        t1ro=-t1*roa_rob*(1.d0+dc*con23*rom13)/t1den


        omega=dexp(-cc*rom13)/(denom*ro**con113)
        omegaro=con13*cc*omega*rom43-
     &    omega*(con113*ro13+con103*dc)/(denom*ro43)

        delta=rom13*(cc+dc/denom)
        deltaro=(-con13*rom43)*(cc+dc/denom)+
     &    (rom13*dc*dc*con13*rom43)/(denom*denom)

c==============================================================
        fcomm=(roa*gro2a+rob*gro2b)*roinv/9.0d0
        fcon=(2.d0**con113)*cf
        fq1=(47.0d0-7.0d0*delta)*gro/18.0d0
        fq2=(2.50d0-delta/18.0d0)
        fq3=(delta-11.0d0)
        f=fcon*(roa**con83+rob**con83)+
     &    fq1*gro-
     &    fq2*(gro2a+gro2b)-
     &    fq3*fcomm  

        fro=(-7.0d0/18.0d0)*deltaro*gro2+
     &    (1.0d0/18.0d0)*deltaro*(gro2a+gro2b)-
     &    deltaro*fcomm+fq3*fcomm*roinv 

        frocomp1=fcon*con83
        frocomp2=-fq3*roinv/9.0d0

        fgro=2.d0*fq1 

        fgrocomp1=-2.d0*fq2
        fgrocomp2=-2.d0*fq3*roinv/9.0d0
c==============================================================
        hcomm=con23*ro*ro
        h=-hcomm*gro2+
     &         (hcomm-roa*roa)*gro2b+(hcomm-rob*rob)*gro2a

        hcommro=con43*ro
        hro=-hcommro*(gro2-gro2b-gro2a)

        hgro=-2.d0*hcomm*gro
c==============================================================

         ab_omega=ac*bc*omega

         comm=t1ro-ac*bc*omegaro*(roa_rob*f+h)-
     &     ab_omega*(roa_rob*fro+hro)

         pota=comm+t1*rob-ab_omega*(rob*f
     &        +roa_rob*(frocomp1*roa**con53+
     &                   frocomp2*gro2a)-2.d0*roa*gro2b)

         potb=comm+t1*roa-ab_omega*(roa*f
     &        +roa_rob*(frocomp1*rob**con53+
     &                   frocomp2*gro2b)-2.d0*rob*gro2a)

         grcompa=-ab_omega*(roa_rob*fgro+hgro) 

         gracompa=-ab_omega*(roa_rob*(fgrocomp1*groa+
     &                                 fgrocomp2*roa*groa)+
     &         2.d0*groa*(hcomm-rob*rob))

         grbcompb=-ab_omega*(roa_rob*(fgrocomp1*grob+
     &                                 fgrocomp2*rob*grob)+
     &         2.d0*grob*(hcomm-roa*roa))

         grcompb=grcompa 
      else
         pota=0.d0
         potb=0.d0
         grcompa=0.d0
         grcompb=0.d0
         gracompa=0.d0
         grbcompb=0.d0
      end if

      return
      end
