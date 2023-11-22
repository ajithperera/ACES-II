#include "flags.h"
      Subroutine Do_FullNR()

      Logical Converged,Bca,Bfc,Occnums_file,Move
      Character*80 Fname 
      Integer Icycle 
      Logical porbrot_in
      Parameter(Max_cycle=20)
#include "flags.com"
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)"), "@-do_FullNR  Entering pCCD full NR  block"
#endif
      
      Converged = .FALSE.
      Icycle = 1
      ImacroItr=1
      ImicroItr=1
!      call putrec(20,"JOBARC","MACRONUM",1,Icycle)
!      call putrec(20,"JOBARC","MICRONUM",1,Icycle)
!      CALL putrec(20,"JOBARC","MACCONVG",1,.False.)
!      CALL putrec(20,"JOBARC","MICCONVG",1,.False.)
!      call putrec(20,"JOBARC","GLOBNORM",1,0.0d0)
!      Do while (.NOT. Converged)
        do jj=1,Max_cycle
          If (Icycle .EQ. 1) Then
             Call Runit("xvmol")
             Call Runit("xvmol2ja")
             Call Runit("xvscf")
          Else
             Call Runit("xvscf")
          Endif
          Call Runit("xvtran")
          Call Runit("xintprc")
          Call Runit("xpccd")

          Icycle = Icycle + 1
          Write(6,"(a,i2,2a)") "@-Do_FullNR starting ",Icycle, " NR",
     +                         " pCCD iteration."
          IF ( Icycle .Gt. Max_cycle)
     +        Write(6,"(a,a,i2,a)") "@-Do_FullNR The maximum allowed",
     +                              " NR pCCD iterations ", Max_cycle,
     +                              " has reached."
      Enddo

      Return
      End
