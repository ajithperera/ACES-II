#include "flags.h"

      program ringcl
      implicit none

      integer iuhf
      logical solve_4lambda,Density 

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
      Call aces_init(icore, i0,icrsiz, iuhf, .true.)
C
      Solve_4lambda = .True.
      Density = (IFLAGS(h_IFLAGS_density) .EQ. 1) 
      If (Density) Solve_4lambda = .True. 

C Respnse density always need left hand wavefunction. If only 
C excited states are requested then turn off lambda. 

      If((IFLAGS(h_IFLAGS_excite) .EQ.  1 .OR. 
     +   IFLAGS(h_IFLAGS_excite)  .EQ.  2 .OR. 
     +   IFLAGS(h_IFLAGS_excite)  .EQ.  4 .OR.
     +   IFLAGS(h_IFLAGS_excite)  .EQ. 13 .OR.
     +   IFLAGS(h_IFLAGS_excite)  .EQ. 14 .OR.
     +   IFLAGS(h_IFLAGS_excite)  .EQ. 15 .OR.
     +   IFLAGS(h_IFLAGS_excite)  .EQ. 16 .OR.
     +   IFLAGS(h_IFLAGS_excite)  .EQ. 17 .OR.
     +   IFLAGS(h_IFLAGS_excite)  .EQ. 18 .OR.
     +   IFLAGS(h_IFLAGS_excite)  .EQ. 19 .OR.
     +   IFLAGS(h_IFLAGS_excite)  .EQ. 20).AND. 
     +   (.NOT. Density)) Solve_4lambda = .False.

      Call Rcl_hbar(Icore(i0),icrsiz/iintfp,Iuhf,Solve_4lambda) 

      Call rcl_driver(Icore(i0),icrsiz/iintfp,Iuhf,Solve_4lambda)

#ifdef _DEBUG_LVL0
CSSS        Call rcl_check_rcchbar(Icore(i0),icrsiz/iintfp,Iuhf) 
#endif

      Call aces_fin
C
c ----------------------------------------------------------------------
      Stop
      End

