      Subroutine Freq_anlytcl(Calc_level)
c
#include "aces2.h"
c
      Character*4 Calc_level
     
      If (Calc_level .EQ. "SCF " .and. .not. ks_scf) Then
         Call Scf_Hess
      Else
         Print*, "The analytic Hessian for post-SCF is not available"
         Call Errex
      Endif
c
      Return
      End
c
      Subroutine Scf_Hess
c
#include "aces2.h"
c     
      Call Scf_eneg
      Call Runit('xvtran')
      Call Runit('xintprc')
      Call Scf_anlytcl_grad(.True.) 
      Call Runit('xcphf')
CSSS      Call Runit('rm -r OPTARC') 
      Call Runit('xjoda')
c  
      Return
      End
