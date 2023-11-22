      Subroutine Post_scf_higher_order_prps(DO_NMR_SHIFTS,
     &                                      DO_NMR_SPNSPN,
     &                                      DO_NLO_PROPS)
c
       Logical DO_NMR_SHIFTS, DO_NMR_SPNSPN, DO_NLO_PROPS 
c
#include "flags.h" 
#include "aces2.h"
c
      IF (DO_NMR_SHIFTS) Then
c
         Call Runit('xlambda')
         Call Runit('xdens')
         Call Runit('xvdint')
         Call Runit('xcphf')
         Call Runit('xnmr')
c
        If (iflags(h_IFLAGS_treat_pert).ne.0) then
c
c o sequential treatment of B-field components requires
c   reevaluation of GIAO integrals
c
            call runit('rm -f I*X')
            call runit('xvdint')
            call runit('xnmr')
            call runit('rm -f I*Y')
            call runit('xvdint')
            call runit('xnmr')
            call runit('rm -f I*Z')
         End if
         Call Runit('xjoda')
c
      Else if (DO_NMR_SPNSPN .OR. DO_NLO_PROPS) Then
c
         Call Runit('xlambda')
         Call Runit('xvprops')
         Call Runit('xvcceh')
c
      Endif
c
      Return
      End
