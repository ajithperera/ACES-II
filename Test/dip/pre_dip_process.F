#include "flags.h"
      Program pre_dip_process
       
      ImpliciT Double Precision (A-H, O-Z)

      Logical act_spc_cc
      Character*4 act
      Dimension Nocc(16)

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "syminf.com"
#include "files.com"

      Luout = 6
C
      call aces_init(icore,i0,icrsiz,iuhf,.True.)
      act_spc_cc = .true.
      act = "READ"
      if (act_spc_cc .and. act .eq."READ")
     &    call modf_amps(ICORE(I0),Icrsiz,IUHF,0,.TRUE.,'T',act,"AMPS")

C
c Read in the new occ numbers or other orbital manipulations from the
C GUESS file. 
c
      Write(6,*) "I am Here"
      Call Readgs(Nocc, Nirrep, Iuhf)
      Call Putrec(20,'JOBARC','OCCUPYA0',Nirrep,NOCC(1))
      If (Iuhf .gt. 0) Call Putrec(20,'JOBARC','OCCUPYB0',Nirrep,
     &                             NOCC(9))
C
C Set the scf_maxcyc to zero and rerun the SCF to built n+2 Fock matrix
C
      iflags(h_IFLAGS_scf_maxcyc) = 0
      Call Putrec(20, "JOBARC", "IFLAGS  ", 100, Iflags)
      Call Runit("xvscf") 
C

      Call aces_fin
       
      End
