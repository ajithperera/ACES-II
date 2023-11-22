
#ifndef _ENERG_COM_
#define _ENERG_COM_

c Analytic values determined in the anlytc subroutine

c All HF values are determined with the HF density
c All correlated values are determined using the correlated relaxed density

c ehfkin   : Hartree-Fock kinetic energy
c eckin    : correlated
c ehfnatr  : Hartree-Fock nuclear-electron attraction energy
c ecnatr   : correlated
c ehfcoul  : Hartree-Fock coulomb energy
c eccoul   : correlated
c ehfx     : Hartree-Fock exchange energy
c ecx      : correlated (using the correlated relaxed density in the
c            Hartree-Fock exchange energy expression)

c ehar     : Hartree energy using HF orbitals
c echar    : Hartree energy using natural orbitals and occupations
c escf     : SCF energy
c ecor     : energy determined using the relaxed density in the HF
c              energy expression

      M_REAL
     &    ehfkin,eckin,ehfnatr,ecnatr,ehfcoul,eccoul,ehfx,ecx,
     &    ehar,echar,escf,ecor
      common /exch/ ehfkin,eckin,ehfnatr,ecnatr,ehfcoul,eccoul,
     &    ehfx,ecx,ehar,echar,escf,ecor
      save /exch/

#endif /* _ENERG_COM_ */

