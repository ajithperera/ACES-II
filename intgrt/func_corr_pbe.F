c     Perdew-Burke-Ernzerhof Correlation Energy Functional
c
c     By S. Ivanov using the original version by Kieron Burke.
c
c     Phys. Rev. Lett. 77, 3865 (1996); ibid 78, 1396 (1997)
c     Phys. Rev. B 54, 16 533 (1996).
c     Phys. Rev. B 45, 13 244 (1992).

      subroutine func_corr_pbe(func)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &    func

      M_REAL
     &   con13,con23,con16,beta,gamma,delta,pi

      parameter (con13=1.d0/3.d0,con23=2.d0/3.d0,con16=1.d0/6.d0)
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(beta=0.06672455060314922d0,delta=beta/gamma)

      M_REAL
     &   confk,consk,
     &   g,g3,sk,twoksg,t,t2,t4,b,b2,ec,fk,pon,
     &   rs2,rs3,q4,q5,h

c g=phi(zeta)
c rs=(3/(4pi*rho))^(1/3)=local Seitz radius
c sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
c twoksg=2*Ks*phi
c t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
c ec=lsd correlation energy
c h=gradient correction to correlation energy

 
      callstack_curr='FUNC_CORR_PBE'  

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

        call func_corr_lsd2(ec)

c PBE correlation energy

        g3=g*g*g
        pon=-ec/(g3*gamma)
        b=delta/(dexp(pon)-1.d0)
        b2=b*b
        t2=t*t
        t4=t2*t2
        rs2=rs*rs
        rs3=rs2*rs
        q4=1.d0+b*t2
        q5=1.d0+b*t2+b2*t4
        h=g3*(beta/delta)*dlog(1.d0+delta*q4*t2/q5)

        func=ro*(ec+h)
      else
        func=0.d0

      end if

      return
      end
