
C JDW 5/12/95 Force calculation in POLYRATE jobs
C     5/25/95 Only let this happen if we have cartesian coordinates.
C             Note that we can in principle do a force with dropped
C             core, we cannot do hessian, hence dropped core addition
C             here is currently of minor value.

      subroutine polyrate

#include "aces2.h"

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

      ndrgeo=0
      call a2getrec(20,'JOBARC','NDROPGEO',1,ndrgeo)

      if (.not.direct) call runit(integral_package)
      call runit('xvmol2ja')
      if (nddo_guess) call runit('xnddo')
      If (dkh .and. contract) Call Runit("xdkh")
      call runit('xvscf')
      if (iflags(2).gt.0) then
         if (ndrgeo.ne.0) then
            call runit('cp JOBARC JOBARC_AM')
            call runit('cp JAINDX JAINDX_AM')
         end if
         call runit('xvtran')
         call runit('xintprc')
         call runit('xvcc')
         call runit('xlambda')
         if (ndrgeo.ne.0) call runchg
         call runit('xdens')
         call runit('xanti')
         call runit('xbcktrn')
      end if
      call runit('xvdint')

      return
      end

