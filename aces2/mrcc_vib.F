
      subroutine mrcc_vib

#include "aces2.h"
#include "parallel_aces.com"

      integer ires1, ires2, itype, istat

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

c YAU - real restarts should be added to mrcc vib freqs

 1000 continue
         if (analytical_gradient) then
            call mrcc_grad
         else
            call energy_4mrcc
         end if
         call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
         if (istat.eq.0) then
            call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
            call runit(xjoda)
         else
#ifdef _HAVE_MPI
            if (nprocs.ne.1) call sync_fd
#endif
            if (irank.eq.0) then
               call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
               call runit(xjoda)
            end if
            return
         end if
      goto 1000

      end

