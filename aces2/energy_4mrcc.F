
      subroutine energy_4mrcc

#include "aces2.h"
#include "flags.h"

      logical brueckner,dkhess           ! ...Watson... !
      integer iter, maxiter, izero, itest, iflags0(100),
     $     iflags2_0(500), iuhf0

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif

      call a2getrec(20,'JOBARC','IFLAGS  ',100,iflags0)
      call a2getrec(20,'JOBARC','IFLAGS2 ',500,iflags2_0)

      if (.not.direct) call runit(integral_package)
      call runit('xvmol2ja')
      iuhf0 = iflags0(11)

      if (mrcc.or.
     &    (iflags(18).ne.0.and.iflags(18).ne.3.or.iflags(23).ne.0.or.
     &     iflags(24).ne.0.or. iflags(25).ne.0.or.iflags(87).ne.0)
     &   ) call runit('xvprops')
C
C
C   ...Watson...
C
      dkhess = (iflags2(h_IFLAGS2_oomp_cyc) .gt. 0)
      if (mrcc .and. (.not. dkhess)) call runit('rm -f VPOUT')
C      if (mrcc) call runit('rm -f VPOUT')
C
C   ...Watson...
C
      if (hf_scf) then
         if (nddo_guess) call runit('xnddo')
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')
      else
         If (dkh .and. contract) Call Runit("xdkh")
         call runit('xvscf')
         call runit('xintgrt')
      endif

      if (fno) then
         call runfno
      end if

      if ( mrcc .or. (iflags(2).gt.0) .or. (iflags(87).ne.0) ) then
         if (dirmp2) then
            call runit('xdirmp2')
         else
            call runit('xvtran')
            call runit('xintprc')
         end if
         if (mrcc) then
            brueckner = iflags(22) .eq. 1
            call runit('xmrcc')
            call a2getrec(20,'JOBARC','IFLAGS  ',100,iflags)
c
            if (iflags(22) .eq. 1) then
c
               if (.not. brueckner) then
c
c brueckner was switched on by mrcc program
c
c     integrals are no longer there
c
                  if (.not.direct) call runit(integral_package)
                  call runit('xvmol2ja')
                  call runit('xvprops')
                  call runit('rm -f VPOUT')
               endif
c
c     Continue brueckner calculation
c
               maxiter = iflags(7)
               izero = 0
               call a2putrec(20, 'JOBARC','BRUKITER', 1, izero)
               call a2putrec(20, 'JOBARC','BRUKTEST', 1, izero)
c
 100           continue
               call a2getrec(20,'JOBARC','BRUKTEST', 1,itest)
               call a2getrec(-1, 'JOBARC', 'BCC_ITER', 1, iter)
c
               if((itest .eq. 0) .and. iter .le. maxiter) then
                  if (nddo_guess) call runit('xnddo')
                  call runit('xvscf')
                  call runit('xvtran')
                  call runit('xintprc')
                  call runit('xmrcc')
                  goto 100
               else
                  write(6,*) ' @aces2: Brueckner calculation all done'
               endif
c
            endif
c
            call rmfiles3b
         endif
      end if
c
c restore original flags
c
      call a2putrec(20,'JOBARC','IFLAGS  ',100,iflags0)
      call a2putrec(20,'JOBARC','IFLAGS2 ',500,iflags2_0)

      return
      end

