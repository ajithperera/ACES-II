      Subroutine Form_new_amps(Tamps,Work,Iwork,Maxcor,Imaxcor,
     +                         Tln,T1ln,T2ln,T2ln_aa,T2ln_bb,T2ln_ab,
     +                         T1ln_aa,T1ln_bb,Iuhf)
 
      Implicit Double Precision(A-H,O-Z)
      Integer T2ln,T2ln_aa,T2ln_bb,T2ln_ab
      Integer Tln,T1ln,T1ln_aa,T1ln_bb
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

      Dimension Tamps(Tln)
      Dimension Work(Maxcor)
      Dimension IWork(Imaxcor)

      List1 = 90
      T1off = Ione + T2ln
      Call Form_new_t1amps(Tamps(T1off),Work,Iwork,Maxcor,Imaxcor,
     +                     LIst1,T1ln,T1ln_aa,T1ln_bb,Iuhf)

      List2 = 63
      Call Form_new_t2amps_aa(Tamps,Work,Iwork,Maxcor,Imaxcor,List2,
     +                        T2ln,T2ln_aa,Iuhf)

      List2 = 66
      Call Form_new_t2amps_ab(Tamps,Work,Iwork,Maxcor,Imaxcor,List2,
     +                        T2ln,T2ln_aa,T2ln_bb,Iuhf)

      Return
      End
