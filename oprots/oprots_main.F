#include "flags.h"

      Program  oprots 
      Implicit none

      Integer iuhf

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
      Call Aces_init(icore, i0,icrsiz, iuhf, .true.)
      If (Iuhf .Eq. 1) Then
         Write(6,*)
         Write(6,"(2a)") " UHF reference orbital rotations are not",
     +                   " supported changed the refrence to RHF"
         Write(6,"(a)")  " and rerun to proceed."
         Call Errex()
      Endif

C Orbital pair rotations
      Call Oprots_driver(Icore(i0),icrsiz/iintfp,Iuhf)

      Call aces_fin
C
c ----------------------------------------------------------------------
      Stop
      End

