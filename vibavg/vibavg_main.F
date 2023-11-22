#include "flags.h"

      program vibavg
      implicit none

      integer iuhf

      External Tdee_oed_pes

c COMMON BLOCKS
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"

c ----------------------------------------------------------------------
C We can not call aces_init directly because aces_int try to initiate
C the chemical system by default and that can fail for finite difference
C frequency calculations (if the last point happens to be in lower 
C symmetry)
CSSS      Call aces_init(icore, i0,icrsiz, iuhf, .true.)

      Call Aces_init_rte
      Call Aces_com_parallel_aces
      Call Aces_ja_init
      Call Getrec(1,'JOBARC','IFLAGS', 100,iflags)
      Call Getrec(1,'JOBARC','IFLAGS2',500,iflags2)

      Icrsiz = Iflags(h_IFLAGS_mem)
      icore(1) = 0
      do while ((icore(1).eq.0).and.(icrsiz.gt.1000000))
         call aces_malloc(icrsiz,icore,i0)
         if (icore(1).eq.0) icrsiz = icrsiz - 1000000
      end do

      Call Vibavg_driver(Icore(i0), icrsiz/iintfp, Iuhf)

      Call aces_fin
C
c ----------------------------------------------------------------------
      stop
      end

