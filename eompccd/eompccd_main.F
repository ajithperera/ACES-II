#include "flags.h"

      Program  eompccd_main
      Implicit none
      Integer iuhf
      Double Precision Scale 

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "sympop.com"
#include "machsp.com"
       
c ----------------------------------------------------------------------
      Call Aces_init(icore, i0,icrsiz, iuhf, .true.)
      If (Iuhf .Eq. 1) Then
         Write(6,*)
         Write(6,"(2a)") " UHF reference pCCD calculations are not",
     +                   " supported changed the refrence to RHF"
         Write(6,"(a)")  " and rerun to proceed."
         Call Errex()
      Endif 

      Write(6,*) 
      Write(6,"(5x,2a)") "------Eom-EE/IP/EA/DIP/DEA pCCD calculations",
     +                   " are performed------"
      Write(6,*) 
      Call Eompccd_hbar(Icore(i0),icrsiz/iintfp,Iuhf)
CSSS      Call Eompccd_ee(Icore(i0),Icrsiz/iintfp,Iuhf) 
      Call Eompccd_ea(Icore(i0),Icrsiz/iintfp,Iuhf) 

      Call aces_fin
c ----------------------------------------------------------------------
      Stop
      End

