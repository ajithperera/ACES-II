#ifndef _CCDLIGHT_VARS.COM
#define _CCDLIGHT_VARS.COM
#include "maxbasfn.par" 

      Integer cc_maxcyc
      Logical Ring_cc
      Double Precision ocn_oa,Ocn_ob,Ocn_va,Ocn_vb
      Double Precision Denom_tol  
      Dimension E_corr(0:500)

      Common /ccdlight_vars/Ring_cc,cc_conv,cc_maxcyc,
     +                      ocn_oa(Maxbasfn),ocn_ob(Maxbasfn),
     +                      ocn_va(Maxbasfn),ocn_vb(Maxbasfn),
     +                      E_corr,Denom_tol

#endif  /* _CCDLIGHT_VARS.COM__ */


