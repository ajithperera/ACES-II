c       LSD correlation potential needed for PBE
c
c       By S. Ivanov, January


      subroutine pot_corr_lsd2(pota,potb,ec,ecrs,eczeta)

      implicit none
#include <aces.h>
#include <dftfunc.com>
#include <callstack.com>

      M_REAL
     &   pota,potb,ec,ecrs,eczeta 

      M_REAL
     &   con13,con43,fzz,gam

       M_REAL
     & rtrs,eu,eurs,ep,eprs,alfm,alfrsm,f,z4,alfc,fz,comm

      parameter (con13=1.d0/3.d0,con43=4.d0/3.d0)

c  In PBE some of the constants are more accurate.
c      parameter(gam=0.5198420997897463295344212145565d0)
c      parameter(fzz=8.d0/(9.d0*gam))

c  In PW91 less accurate values for gam and fzz have been used.
      parameter(gam=0.5198421d0,fzz=1.709921d0)

c eu=unpolarized LSD correlation energy
c eurs=deu/drs
c ep=fully polarized LSD correlation energy
c eprs=dep/drs
c alfm=-spin stiffness
c alfrsm=-dalpha/drs
c f=spin-scalining factor 
c construct ec
c ecrs=dec/drs
c eczeta=dec/dzeta
      
      callstack_curr='POT_CORR_LSD2'

      rtrs=dsqrt(rs)

      call gcor2(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     &    0.49294D0,rtrs,eu,eurs)
      call gcor2(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     &    0.62517D0,rtrs,ep,eprs)
      call gcor2(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     &    0.49671D0,rtrs,alfm,alfrsm)
      alfc=-alfm

      z4=zeta**4
      f=((1.d0+zeta)**con43+(1.d0-zeta)**con43-2.d0)/gam
      ec=eu*(1.d0-f*z4)+ep*f*z4-alfm*f*(1.d0-Z4)/fzz

      ecrs=eurs*(1.d0-f*z4)+eprs*f*z4-alfrsm*f*(1.d0-z4)/fzz
      fz=con43*((1.d0+zeta)**con13-(1.d0-zeta)**con13)/gam
      eczeta=4.d0*(zeta*zeta*zeta)*f*(ep-eu+alfm/fzz)+
     &  fz*(z4*ep-z4*eu-(1.d0-z4)*alfm/fzz)

      comm=ec-rs*ecrs/3.d0-zeta*eczeta

      pota=comm+eczeta
      potb=comm-eczeta

      return
      end
