      Subroutine Pdcc_fixf(Work,Maxcor,Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Maxcor)

#include "parcc.com"

      Call Dcc_hbar_fae(Work,Maxcor,Iuhf)
      Call Dcc_hbar_fmi(Work,Maxcor,Iuhf)
      Call Dcc_fixfbar(Work,Maxcor,Iuhf)

      REturn
      End 
