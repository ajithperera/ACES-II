      program pdens_main
      implicit none
#include "flags.h"
c
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
      integer iuhf
c ----------------------------------------------------------------------
      call aces_init_rte
      call aces_com_parallel_aces
      call aces_ja_init
      call getrec(1,'JOBARC','IFLAGS', 100,iflags)
      call getrec(1,'JOBARC','IFLAGS2',500,iflags2)
c ----------------------------------------------------------------------

      icrsiz = iflags(h_IFLAGS_mem)
      icore(1) = 0
      do while ((icore(1).eq.0).and.(icrsiz.gt.1000000))
         call aces_malloc(icrsiz,icore,i0)
         if (icore(1).eq.0) icrsiz = icrsiz - 1000000
      end do
      if (icore(1).eq.0) then
         print *, '@MAIN: unable to allocate at least ',
     &            1000000,' integers of memory'
         call aces_exit(1)
      end if

c ----------------------------------------------------------------------

      call aces_io_init(icore,i0,icrsiz,.true.)

      if (iflags(h_IFLAGS_ref) .eq. 0) then
         iuhf = 0
      else
        iuhf = 1
      endif
c
#ifdef _DEBUG_LVL0
         print*, "Entering PERT_DENS_MAIN" 
         print*, "The iuhf:", iuhf
#endif
      Write(6,*)
      Write(6,"(T2,a)") "Perturbed densities are computed"

         call aces_init_chemsys
         call pert_dens_main(icore(i0), icrsiz, iuhf)
c ----------------------------------------------------------------------

 9997 continue
      call aces_io_fin
 9998 continue
c      call c_free(icore)
 9999 continue
      call aces_ja_fin
C
      Write(6,*)
      Write(6,"(T2,a)") "Perturbed densities were successfully computed"
      call c_exit(0)
      end

