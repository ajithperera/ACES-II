
      subroutine opt_num

#include "aces2.h"
#include "parallel_aces.com"

c   o generate an initial approximate Hessian. 12/04, Ajith Perera
      if (mopac_guess) call runit('xmopac')

      istat = 0
      do while (istat.eq.0)
         call energy
#ifdef _HAVE_MPI
         call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
         if (nprocs.ne.1.and.istat.ne.0) call sync_fd
#endif
         call c_gtod(is,ius)
         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
     &                        is+1.d-6*ius-dTimeStart,' seconds'
         call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
         call runit(xjoda)
         call a2getrec(1,'JOBARC','JODADONE',1,istat)
      end do

      return
      end

