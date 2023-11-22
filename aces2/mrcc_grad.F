
      subroutine mrcc_grad

#include "aces2.h"

      logical brueckner
      integer iter, maxiter, izero, itest
c
#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

      ndrgeo=0
      call a2getrec(20,'JOBARC','NDROPGEO',1,ndrgeo)
      brueckner = iflags(22) .eq. 1

      if (.not.direct) call runit(integral_package)
      call runit('xvmol2ja')
      call runit('xvprops')
c     
      if (brueckner) then
         write(6,*) ' @mrcc_grad: run Brueckner calculation '
c     
c     Perform Brueckner iteration
c     
         maxiter = iflags(7)
         izero = 0
         call a2putrec(20, 'JOBARC','BRUKITER', 1, izero)
         call a2putrec(20, 'JOBARC','BRUKTEST', 1, izero)
         call a2putrec(20, 'JOBARC','BCC_ITER', 1, izero)
         call a2putrec(20, 'JOBARC','PREVBRUK', 1, izero)
         call a2putrec(20, 'JOBARC','NDIMBDIS', 1, izero)
c     
 100     continue
         call a2getrec(20,'JOBARC','BRUKTEST', 1,itest)
         call a2getrec(-1, 'JOBARC', 'BCC_ITER', 1, iter)
c     
         if((itest .eq. 0) .and. iter .le. maxiter) then
            if (nddo_guess) call runit('xnddo')
            If (dkh .and. contract) Call Runit("xdkh")
            call runit('xvscf')
            if (ndrgeo.ne.0) then
               call runit('cp JOBARC JOBARC_AM')
               call runit('cp JAINDX JAINDX_AM')
            end if
            call runit('xvtran')
            call runit('xintprc')
            call runit('xmrcc')
            goto 100
         else
            write(6,*) ' @aces2: Brueckner calculation all done'
         endif
         call a2getrec(20,'JOBARC','BRUKTEST', 1,itest)
         if (itest .ne. 1) then
            write(6,*) ' Brueckner calculation did not converge '
            write(6,*) ' Abort run '
            call errex
         endif
      else
         if (nddo_guess) call runit('xnddo')
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')
         if (ndrgeo.ne.0) then
            call runit('cp JOBARC JOBARC_AM')
            call runit('cp JAINDX JAINDX_AM')
         end if
c     
         call runit('xvtran')
         call runit('xintprc')
         call runit('xmrcc')
      endif                     ! brueckner branch
c     
      if (ndrgeo.ne.0) call runchg
      call runit('xdens')
      call runit('xanti')
      call runit('xbcktrn')
      call runit('xvdint')
      if (iflags2(3) .eq. 1) then
c     mn Resonance raman calculation
         call runit('xjoda')
      endif
      call rmfiles3
      return
      end

