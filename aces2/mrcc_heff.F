
      subroutine mrcc_heff

#include "aces2.h"
#include "flags.h"
#include "parallel_aces.com"

      integer istat

      write(*,*) '@ACES2: Entered ', "mrcc_heff.F, irank", irank

      istat=0
      do while (istat.eq.0)
         call energy_4mrcc
         write(6,*) ' Done energy_4mrcc, back to xaces2'
         call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
              if (istat.eq.0) then
                 call c_gtod(is,ius)
                 print '(a,f10.1,a)',
     &                 'ACES2: Total elapsed time is ',
     &                     is+1.d-6*ius-dTimeStart,' seconds'
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
             end if
                call rmfiles
      end do
      
      Return
      End 


