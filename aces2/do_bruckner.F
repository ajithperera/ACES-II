#include "flags.h"
      Subroutine do_brueckner

      Logical Converged
      Logical Rot_grad
      Integer Icycle 
C
#include "aces2.h"

#ifdef _DEBUG_LVL0
      Write(6,"(a)"), "@-DO_POST_Brueckner Entering Brueckner block"
#endif
      Rot_grad = Iflags2(h_IFLAGS2_oo_method) .Eq. 4
      Converged = .FALSE.
      Icycle = 1
      Do while (.NOT. Converged)
          If (Icycle .EQ. 1) Then
              Call Scf_eneg
          Else
              If (dkh .and. contract) Call Runit("xdkh")
              Call Runit("xvscf")
          Endif
          Call Change_orbital_space
          Call Prep4_post_scf
          Call Post_scf_eneg_4brueck(Rot_grad)
          Call A2getrec(20, 'JOBARC', 'BRUKTEST', 1, Itest)
          IF (Itest .EQ. 1) Converged = .TRUE.
          Icycle = Icycle + 1
      Enddo
c
      if (First_order_props) Then
          Call Runit("xlambda")     
          Call Runit("xdens")     
          Call Runit("xvprops")     
      Endif
C
C Excited state calculations with Brueckner orbitals. At convergence
C we have T2 amplitudes to form Hbar.
C
      If (Iflags(h_IFLAGS_excite) .GT. 2 .AND.
     &         Iflags2(h_IFLAGS2_copt) .EQ. 1) Then
         Call Runit("xhbar")
         Call Runit("xfsip")
      Else If (Iflags(h_IFLAGS_excite) .GT. 2 .AND. (.NOT.
     &         Iflags2(h_IFLAGS2_copt) .EQ. 1)) Then
         Call Runit("xhbar")
         Call Runit("xvee")
      Endif 
      
      Return
      End
