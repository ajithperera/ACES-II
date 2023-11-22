c     Perdew-Wang 91 Correlation Potential
c
c     By S. Ivanov using the original version by Kieron Burke.
c
c     Phys. Rev. B 54, 16 533 (1996).
c     Phys. Rev. B 45, 13 244 (1992).


      subroutine pot_corr_pw91(pota,grcompa,potb,grcompb)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    pota,potb,grcompa,grcompb
 
      M_REAL
     &   con13,con23,con16,pi

      parameter (con13=1.d0/3.d0,con23=2.d0/3.d0,con16=1.d0/6.d0)

      M_REAL
     &   xnu,cc0,cx,
     &   c1,c2,c3,c4,c5,c6,a4,
     &   alfa,beta,delta,eta

      parameter (xnu=15.75592d0,cc0=0.004235d0,cx=-0.001667212d0)
      parameter (c1=0.002568d0,c2=0.023266d0,c3=7.389d-6,c4=8.723d0,
     &           c5=0.472d0,c6=7.389d-2,a4=100.0d0)
      parameter(alfa=0.090d0,beta=xnu*cc0,delta=2.d0*alfa/beta)
      parameter (eta=1.d-12)

      M_REAL
     &   confk,consk,
     &   g,g3,g4,sk,twoksg,t,t2,t4,t6,b,b2,ec,fk,pon,
     &   rs2,rs3,q4,q5,q6,q7,r0,r1,r2,r3,r4,cc,coeff,h0,h1,h

      M_REAL
     &   ecrs,eczeta,gz,ccrs,rs13,bg,bec,fac,q8,q9,
     &   h0b,h0z,h0rs,h0t,h1z,h1rs,h1t,hz,hrs,ht,comm

c g=phi(zeta)
c rs=(3/(4pi*rho))^(1/3)=local Seitz radius
c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
c twoksg=2*Ks*phi
c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
c ec=lsd correlation energy
c h=gradient correction to correlation energy

      callstack_curr='POT_CORR_PW91'  

      pi = acos(-1.d0)

      if (ro.ge.thresh) then
        
        confk=(3.d0*pi*pi)**con13
        consk=dsqrt(4.d0*confk/pi)

        g=((1.d0+zeta)**con23+(1.d0-zeta)**con23)*0.5d0
        fk=confk*(ro**con13)
        sk=consk*(ro**con16)
        twoksg=2.d0*sk*g
        t=gro/(twoksg*ro)
     

c find LSD energy contributions

        call pot_corr_lsd2(pota,potb,ec,ecrs,eczeta)

c PW91 correlation potential

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

        t6=t4*t2

        gz=(((1.d0+zeta)**2+eta)**(-con16)-
     &   ((1.d0-zeta)**2+eta)**(-con16))/3.d0

        ccrs=(c2+2.d0*c3*rs)/q7-q6*(c4+2.d0*c5*rs+3.d0*c6*rs2)/q7**2
        rs13=rs*con13
        r4 =rs13*ccrs/coeff
        fac=delta/b+1.d0
        bg=-3.d0*b2*ec*fac/(beta*g4)
        bec=b2*fac/(beta*g3)
        q8=q5*q5+delta*q4*q5*t2
        q9=1.d0+2.d0*b*t2

        h0b=-beta*g3*b*t6*(2.d0+b*t2)/q8
        h0rs=-rs13*ecrs*h0b*bec

        h0z=3.d0*gz*h0/g+h0b*(bg*gz+bec*eczeta)
        h0t=2.d0*beta*g3*q9/q8

        h1rs=r3*r2*t2*(-r4+r1*t2*con13)
        h1z=gz*r3*r2*t2*(3.d0-4.d0*r1*t2)/g
        h1t=2.d0*r3*r2*(1.d0-r1*t2)

        hrs=h0rs+h1rs
        hz=h0z+h1z
        ht=h0t+h1t


        comm=h+hrs-hz*zeta
     &    +(t2*ht*gz*zeta/g)-7.d0*con16*t2*ht

        pota=pota+comm+hz-(t2*ht*gz/g)
        potb=potb+comm-hz+(t2*ht*gz/g)

        comm=ro*ht*(t2/gro)

        grcompa=comm
        grcompb=comm
      else
        pota=0.d0
        potb=0.d0
        grcompa=0.d0
        grcompb=0.d0
      end if

      return
      end
