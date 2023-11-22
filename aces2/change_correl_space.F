      Subroutine Change_correl_space

#include "flags.h"
#include "aces2.h"
c
c KJW 6/16/98 FNO reductions of virtual space. Any other scheme
c that change the occupied/virtual space dimensions can go
c in here.
c
      if (fno) then
        call runfno
        call a2getrec(20,'JOBARC','IFLAGS  ',100,iflags)
        call a2getrec(20,'JOBARC','IFLAGS2 ',500,iflags2)
      end if
c 
      Return
      End
