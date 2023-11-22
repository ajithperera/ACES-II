c     Perdew-Burke-Ernzerhof Correlation Potential
c
c     By S. Ivanov using the original version by Kieron Burke.
c
c     Phys. Rev. Lett. 77, 3865 (1996); ibid 78, 1396 (1997)
c     Phys. Rev. B 54, 16 533 (1996).
c     Phys. Rev. B 45, 13 244 (1992).

      subroutine pot_corr_pbe(pota,grcompa,potb,grcompb)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    pota,potb,grcompa,grcompb

      M_REAL
     &   con13,con23,con16,beta,gamma,delta,eta,pi

      parameter (con13=1.d0/3.d0,con23=2.d0/3.d0,con16=1.d0/6.d0)
      parameter (gamma=0.03109069086965489503494086371273d0)
      parameter (beta=0.06672455060314922d0,delta=beta/gamma)
      parameter (eta=1.d-12)

      M_REAL
     &   confk,consk,
     &   g,g3,sk,twoksg,t,t2,t4,b,b2,ec,fk,pon,
     &   q4,q5,h,ecrs,eczeta

      M_REAL
     &   g4,t6,gz,fac,bg,bec,q8,q9,hb,hrs,qln,
     &   hz,ht,comm

c g=phi(zeta)
c rs=(3/(4pi*rho))^(1/3)=local Seitz radius
c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
c twoksg=2*Ks*phi
c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
c ec=lsd correlation energy
c h=gradient correction to correlation energy

 
      callstack_curr='POT_CORR_PBE'  

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

c PBE correlation potential

        g3=g*g*g
        pon=-ec/(g3*gamma)
        b=delta/(dexp(pon)-1.d0)
        b2=b*b
        t2=t*t
        t4=t2*t2
        q4=1.d0+b*t2
        q5=1.d0+b*t2+b2*t4
        qln=1.d0+delta*q4*t2/q5
        h=g3*(beta/delta)*dlog(qln)

        g4=g3*g
        t6=t4*t2
        gz=(((1.d0+zeta)**2+eta)**(-con16)-
     &    ((1.d0-zeta)**2+eta)**(-con16))*con13

        fac=delta/b+1.d0 
        bg=-3.d0*b2*ec*fac/(beta*g4) 
        bec=b2*fac/(beta*g3) 
        q8=q5*q5+delta*q4*q5*t2 
        q9=1.d0+2.d0*b*t2 
        hb=-beta*g3*B*t6*(2.d0+b*t2)/q8

        hrs=hb*bec*ecrs

        hz=3.d0*gz*h/g+hb*(bg*gz+bec*eczeta)

        ht=2.d0*beta*g3*q9/q8

        comm=h-con13*rs*hrs-hz*zeta
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
