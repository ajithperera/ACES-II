#include "flags.h"
      Subroutine do_brueckner(BCa,Bfc,Occnums_file)

      Logical Converged,Bca,Bfc,Occnums_file,Move
      Character*80 Fname 
      Integer Icycle 
      Parameter(Max_cycle=20)
C
#ifdef _DEBUG_LVL0
      Write(6,"(a)"), "@-do_brueckner Entering Brueckner block"
#endif
      Converged = .FALSE.
      Icycle = 1
      Do while (.NOT. Converged)
          If (Icycle .EQ. 1) Then
             Call Runit("xvmol")
             Call Runit("xvmol2ja")
             Move = (Occnums_file .and. Bfc .and. .not. Bca)
             If (Move) Call Runit("cp OCCNUMS_DUP OCCNUMS")
             Call Runit("xvscf")
          Else
             Call Gfname("OCCNUMS",Fname,Ilength)
             Inquire(File=Fname(1:7),Exist=Occnums_file)
             If (Occnums_file)  Call Runit("rm OCCNUMS")
             Move = (Occnums_file .and. Bfc .and. .not. Bca)
             If (Move) Call Runit("cp OCCNUMS_DUP OCCNUMS")
             Call Runit("xvscf")
          Endif
          If (Occnums_file) Call Runit("cp OCCNUMS_DUP OCCNUMS")
          Call Runit("xccsd_light") 
          Call A2getrec(20, 'JOBARC', 'BRUKTEST', 1, Itest)
          Icycle = Icycle + 1
          Write(6,"(a,i2,2a)") "  Starting ",Icycle, " Brueckner",
     +                         " iteration."
          IF (Itest .EQ. 1 .or. Icycle .Gt. Max_cycle ) 
     +        Converged = .TRUE.
          IF (itest .NE. 1 .AND.  Icycle .Gt. Max_cycle)
     +        Write(6,"(a,a,i2,a)") "  The maximum allowed Brueckner",
     +                              " iterations ", Max_cycle, 
     +                              " has reached." 
      Enddo
      
      Return
      End
