#include "flags.h"
      Subroutine Tdcc_form_mu_dot_s(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"

C First three terms in Eqn. 25 of NDeP paper. The last there terms
C are zero for one-particle perturbation operators.

      Call Tdcc_form_mu_dot_sa(Work,Memleft,Irrepx,Iuhf)

#ifdef _DEBUG_LVL0
      Call tdcc_sa1_r_debug(Work,Memleft,Irrepx,Iuhf)
      Call tdcc_sa2_r_debug(Work,Memleft,Irrepx,Iuhf)
#endif
 
      Return
      End
     
           
        
