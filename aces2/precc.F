
#include "flags.h"

      subroutine precc(ndrgeo)
      implicit none

      integer ndrgeo, ndropa
      integer itest

#include "aces2.h"

c   o vtran/intprc were already run if HFSTAB=(ON|FOLLOW)
c     we have to re-run them if (HFSTAB=FOLLOW && 'SCFKICK'!=0) || DROPMO!=0
      if (iflags(h_IFLAGS_hfstab).eq.2) then
         call a2getrec(-1,'JOBARC','SCFKICK',1,itest)
      else
         itest=0
      end if

c   o use NUMDROPA since ndrgeo is 0 for single-point calcs
      call a2getrec(-1,'JOBARC','NUMDROPA',1,ndropa)

      if (iflags(h_IFLAGS_hfstab).eq.0.or.ndropa.ne.0.or.
     &    (iflags(h_IFLAGS_hfstab).eq.2.and.itest.ne.0)  ) then
         if (ndrgeo.ne.0)  then
            call runit('cp JOBARC JOBARC_AM')
            call runit('cp JAINDX JAINDX_AM')
         end if
         call runit('xvtran')
         if (iflags2(h_IFLAGS2_kucharski).eq.1) then
c         o KJW 9/29/99 for Stan Kucharski's ccsdtq program
            call runit('cp HF2 HF2N')
            call runit('xintprc')
            call runit('mv HF2N HF2')
         else
            call runit('xintprc')
         end if
      end if

      return
      end

