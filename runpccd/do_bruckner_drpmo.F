#include "flags.h"
      Subroutine do_brueckner_drpmo(Xpccd)

      Implicit Double Precision(A-H,O-Z)
      Logical Converged,Bca,Bfc,Occnums_file,Move
      Character*80 Fname 
      Integer Icycle 
      Logical porbrot_in
      Logical Xpccd 
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
C
#ifdef _DEBUG_LVL0
      Write(6,"(2a)"), " @-do_brueckner Entering pCCD Brueckner/or",
     +                 " orbital optimization block"
#endif
      
      Converged = .FALSE.
      Max_cycle = Iflags2(h_IFLAGS2_oo_maxcyc)
      No_diis   = Iflags2(h_IFLAGS2_oo_diis)

      Inquire(file="porbrot.info",exist=porbrot_in)
      Icycle = 1
      Do while (.NOT. Converged)
          If (Icycle .EQ. 1) Then
             Call Runit("xvmol")
             Call Runit("xvmol2ja")
             Call Runit("xvscf")
             If (porbrot_in) Then
                Write(6,"(3a)")" Warning! - porbrot_in file is",
     +                         " present and specified extrenal",
     +                         " orbital rotations are performed!"
                Call Runit("xoprots")
             Endif 
          Else
             Call Runit("xvscf")
          Endif

CSSS          Call Aces_ja_init 
CSSS          Call putrec(20,'JOBARC','IFLAGS  ',100,iflags)
CSSS          Call Aces_ja_fin

          Call prep4_post_scf() 
          Call Runit("xvcc")
          Call Runit("xlambda")
          Call runchg
          Call Runit("xdens")
          If (Xpccd) Then
             Call Runit("xpccd_drpmo")
CSSS             Call Runit("xpsi4dbg")
          Else
             Call Runit("xvcc")
          Endif
          If (Icycle .Eq. 1) Then
             Call A2getrec(20, 'JOBARC', 'SCFENEG ',Iintfp, Eref)
             Call A2getrec(20, 'JOBARC', 'TOTENERG',Iintfp, Ecor)
             Call A2putrec(20, 'JOBARC', 'SCFENEG0',Iintfp, Eref)
             Call A2putrec(20, 'JOBARC', 'TOTENEG0',Iintfp, Ecor) 
             Call Aces_ja_fin
          Endif 

          Call A2getrec(0, 'JOBARC', 'BRUKTEST', Ibtest, Ijunk)
          Call A2getrec(0, 'JOBARC', 'GRADTEST', Igtest, Ijunk)
          If (Ibtest .Gt. 0) Call A2getrec(20, 'JOBARC', 'BRUKTEST',
     +                                    1, Ibtest)
          If (Igtest .Gt. 0) Call A2getrec(20, 'JOBARC', 'GRADTEST',
     +                                    1, Igtest)
          Icycle = Icycle + 1
          If (Ibtest .Eq. 0) Then
          Write(6,"(a,i4,2a)") "  Starting ",Icycle, " Brueckner",
     +                         " pCCD iteration."
          Elseif (Igtest .Eq. 0) Then
          Write(6,"(a,i4,2a)") "  Starting ",Icycle," orbital ",
     +                         "optimization pCCD iteration."
          Endif
C--------- Eliminate after debugging is completed-----------
          Call Aces_ja_init 
          Call putrec(20,'JOBARC','ORBOPITR',1,Icycle)
          Call Aces_ja_fin
C----------------------------------------------------------

          IF (Ibtest .EQ. 1 .or. Igtest .EQ. 1)  Then
             Converged = .TRUE.
          Else if (Icycle .Gt. Max_cycle) Then
          If (Ibtest .Eq. 0) Then
             If (Max_cycle .Eq. 1) Stop
             Write(6,"(a,a,i3,a)") "  The maximum allowed Brueckner",
     +                             " pCCD iterations ", Max_cycle, 
     +                             " has reached and no convergence!."
            Write(6,*)
            Call Errex 
          Elseif (Igtest .Eq. 0) Then
             If (Max_cycle .Eq. 1) Stop
             Write(6,"(a,a,i3,a)") "  The maximum allowed orbital",
     +                             " optimization pCCD iterations ", 
     +                                Max_cycle, 
     +                             " has reached and no convergence!" 
            Write(6,*)
            Call Errex 
          Endif 
          Endif 
      Enddo
      
      Return
      End
