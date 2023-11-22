      Subroutine Form_new_pamps(Tx_peom,Work,Maxcor,Iwork,Imaxcor,
     +                          T1ln,T1ln_aa,T1ln_bb,T2ln,T2ln_aa,
     +                          T2ln_bb,T2ln_ab,Nsizec,Ioffset,
     +                          Irrepx,Iuhf)
 
      Implicit Double Precision(A-H,O-Z)
      Integer T2ln,T2ln_aa,T2ln_bb,T2ln_ab
      Integer T1ln,T1ln_aa,T1ln_bb
      Integer T1off

#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "symoff.com"
#include "active_space.com"

      Data Ione /1/

      Dimension Tx_peom(Nsizec)
      Dimension Work(Maxcor)
      Dimension IWork(Imaxcor)

      List1 = 493
      Ioff  = Ione
      Call Form_new_t1pamps(Tx_peom(Ioff),Work,Iwork,Maxcor,Imaxcor,
     +                      List1,T1ln,T1ln_aa,Ioffset,Irrepx,Iuhf)

      T1off = T1ln_aa + Iuhf*T1ln_bb 
      Ioff  = T1off + Ione
      List2 = 453
      Call Form_new_t2pamps_aa(Tx_peom(Ioff),Work,Iwork,Maxcor,Imaxcor,
     +                         List2,T2ln,T2ln_aa,Ioffset,Irrepx,Iuhf)

      List2 = 456
      Call Form_new_t2pamps_ab(Tx_peom(Ioff),Work,Iwork,Maxcor,Imaxcor,
     +                         List2,T2ln,T2ln_aa,T2ln_bb,Ioffset,
     +                         Irrepx,Iuhf)

      Return
      End
