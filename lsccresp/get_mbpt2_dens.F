      Subroutine Get_mbpt2_dens(T2amps,Doo,Dvv,Dov,Work,Maxcor,T2ln,
     &                          T2ln_aa,T2ln_bb,T2ln_ab,Iuhf)

      Implicit Double Precision (A-H,O-Z)
      Integer T2ln,T2ln_aa,T2ln_bb,T2ln_ab

#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "info.com"

      Dimension T2amps(T2ln)
      Dimension Work(Maxcor)
      Dimension Doo(Nfmi(1)+Iuhf*Nfmi(2))
      Dimension Dvv(Nfea(1)+Iuhf*Nfea(2))
      Dimension DOv(Nt(1)+Iuhf*Nt(2))

      Call Dzero(Doo,Nfmi(1)+Iuhf*Nfmi(2))
      Call Dzero(Dvv,Nfea(1)+Iuhf*Nfea(2))

      Call Densoo_mp2(Doo,T2amps,Work,Maxcor,T2ln,T2ln_aa,T2ln_bb,
     +                T2ln_ab,Iuhf)

      Call Densvv_mp2(Dvv,T2amps,Work,Maxcor,T2ln,T2ln_aa,T2ln_bb,
     +                T2ln_ab,Iuhf)

      Return
      End 
