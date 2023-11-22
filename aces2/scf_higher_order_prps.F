      Subroutine Scf_higher_order_prps(DO_TDHF, DO_NMR_SHIFTS, 
     &                                 DO_NMR_SPNSPN, DO_NLO_PROPS)
c  
      Logical DO_TDHF, DO_NMR_SHIFTS, DO_NMR_SPNSPN,
     &        DO_NLO_PROPS
c
#include "aces2.h"
c
      If (DO_TDHF) Then
         Call Runit('xvprops')
         Call Runit('xvtran')
         Call Runit('xtdhf')
      Else if (DO_NMR_SHIFTS .OR. DO_NMR_SPNSPN) Then
         Call Runit('xvtran')
         Call Runit('xintprc')
         Call Runit('xvdint')
         Call Runit('xcphf')
         Call Runit('xjoda')
      Else if (DO_NLO_PROPS) Then
         Call Runit('xvtran')
         Call Runit('xintprc')
         Call Runit('xvdint')
         Call Runit('xcphf')
      Endif
c
      Return
      End
