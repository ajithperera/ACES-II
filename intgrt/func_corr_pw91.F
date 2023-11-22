c     Perdew-Wang Correlation Energy Functional
c
c     By S. Ivanov using the original version by Kieron Burke.
c
c     Phys. Rev. B 54, 16 533 (1996).
c     Phys. Rev. B 45, 13 244 (1992).

      subroutine func_corr_pw91(func)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &   con13,con23,con16,pi

      parameter (con13=1.d0/3.d0,con23=2.d0/3.d0,con16=1.d0/6.d0)

      M_REAL
     &   xnu,cc0,cx,
     &   c1,c2,c3,c4,c5,c6,a4,
     &   alfa,beta,delta

      parameter (xnu=15.75592d0,cc0=0.004235d0,cx=-0.001667212d0)
      parameter (c1=0.002568d0,c2=0.023266d0,c3=7.389d-6,c4=8.723d0,
     &           c5=0.472d0,c6=7.389d-2,a4=100.0d0) 
      parameter(alfa=0.090d0,beta=xnu*cc0,delta=2.d0*alfa/beta)

      M_REAL
     &   confk,consk,
     &   g,g3,g4,sk,twoksg,t,t2,t4,b,b2,ec,fk,pon,
     &   rs2,rs3,q4,q5,q6,q7,r0,r1,r2,r3,cc,coeff,h0,h1,h


c g=phi(zeta)
c rs=(3/(4pi*rho))^(1/3)=local Seitz radius
c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
c twoksg=2*Ks*phi
c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
c ec=lsd correlation energy
c h=gradient correction to correlation energy

 
      callstack_curr='FUNC_CORR_PW91'  

      pi = acos(-1.d0)

      if (ro.ge.thresh) then
        
        confk=(3.d0*pi*pi)**con13
        consk=dsqrt(4.d0*confk/pi)

        g=((1.d0+zeta)**con23+(1.d0-zeta)**con23)*0.5d0
        fk=confk*(ro**con13)
        sk=consk*(ro**con16)
        twoksg=2.d0*sk*g
        t=gro/(twoksg*ro)
     

c Find the LSD energy contributions
c using more accurate values for some of the constants in
c the LSD. See func_corr_lsd2(ec).

        call func_corr_lsd2(ec)

c PW91 correlation energy

        g3=g*g*g
        g4=g3*g
        pon=-ec*delta/(g3*beta)

        b=delta/(dexp(pon)-1.d0)
        b2=b*b

        t2=t*t
        t4=t2*t2

        rs2=rs*rs
        rs3=rs2*rs

        q4=1.d0+b*t2
        q5=1.d0+b*t2+b2*t4
        q6=c1+c2*rs+c3*rs2
        q7=1.d0+c4*rs+c5*rs2+c6*rs3

        cc=-cx+q6/q7
        coeff=cc-cc0-3.d0*cx/7.d0
 
        r0=0.663436444d0*rs
        r1=a4*r0*g4
        r2=xnu*coeff*g3
        r3=dexp(-r1*t2)

        h0=g3*(beta/delta)*dlog(1.d0+delta*q4*t2/q5)
        h1=r3*r2*t2
        h=h0+h1

        func=ro*(ec+h)
      else
        func=0.d0
      end if

      return
      end
