#include "flags.h"
      Subroutine Tdcc_form_mutilde_dot_s(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"

C Terms in Eqn. 25 of NDeP paper. 

      Call Tdcc_form_mutilde_dot_sa(Work,Memleft,Irrepx,Iuhf)

#ifdef _DEBUG_LVL0
      Call tdcc_sa1_l_debug(Work,Memleft,Irrepx,Iuhf)
      call tdcc_sa2_l_debug(Work,Memleft,Irrepx,Iuhf)
#endif

      Call Tdcc_form_mutilde_dot_sb(Work,Memleft,Irrepx,Iuhf)

#ifdef _DEBUG_LVL0
      Call tdcc_sb_l_debug(Work,Memleft,Irrepx,Iuhf)
#endif

      Call Tdcc_form_mutilde_dot_sc(Work,Memleft,Irrepx,Iuhf)

#ifdef _DEBUG_LVL0
      Call tdcc_sc_l_debug(Work,Memleft,Irrepx,Iuhf)
#endif
 
      Return
      End
     
           
        
