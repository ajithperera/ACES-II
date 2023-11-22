
      subroutine endscf
      implicit none
      integer iRc, iOption

#include "molcas.com"

      call aces_fin

      if (seward) then
c      o $MOLCAS/src/scf/scf.f finalization routines... (fatal)
         iRc = -1
         iOption = 0
         call clsord(iRc,iOption)
         if (iRc.ne.0) then
            write(*,*) '@ENDSCF: There was an error closing the Molcas',
     &                 ' 2-electron integral file.'
            call errex
         end if
cSEWTRACE         call qexit("VSCF")
cSEWARD         call finish(itmp)
      end if

      stop
      end

