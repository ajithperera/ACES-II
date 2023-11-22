
      subroutine opt_numrcl(Calc_level)

#include "aces2.h"
#include "parallel_aces.com"

      Character*4 Calc_level
      Logical Converged
#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

c   o generate an initial approximate Hessian. 12/04, Ajith Perera
      if (mopac_guess) call runit('xmopac')

      istat = 0
      do while (istat.eq.0)
c
         If (Calc_level .Eq. "SCF ") Then
            Call Scf_eneg
         Else
c
c At the moment I do not know whether HFSTABILITY=FOLLOW should be allowed.
c Questions are do we do it only in the first cycle, or every cycle or 
c intermittently.

            If (h_IFLAGS_brueckner .EQ. 1) Then
               Do while (.NOT. Converged)
                  Call Scf_eneg
                  Call Change_orbital_space
                  Call Prep4_post_scf
                  Call Post_scf_eneg 
                  Call A2getrec(20, 'JOBARC', 'BRUKTEST', 1, Itest)
                  IF (Itest .EQ. 1) Converged = .TRUE.
               Enddo
            Else
               Call Scf_eneg
               Call Change_orbital_space
               Call Prep4_post_scf                
               Call Post_scf_eneg
            Endif
         Endif
c
#ifdef _HAVE_MPI
         call a2getrec(1,'JOBARC','LASTGEOM',1,istat)
         if (nprocs.ne.1.and.istat.ne.0) call sync_fd
#endif
CSSS         call c_gtod(is,ius)
CSSS         print '(a,f10.1,a)', 'ACES2: Total elapsed time is ',
CSSS     &                        is+1.d-6*ius-dTimeStart,' seconds'
         call a2putrec(1,'JOBARC','DIRTYFLG',1,0)
         call runit(xjoda)
         call a2getrec(1,'JOBARC','JODADONE',1,istat)
      end do

      return
      end

