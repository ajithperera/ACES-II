
      subroutine energy

#include "aces2.h"

      logical brueckner
      integer iter, maxiter, izero, itest

      if (.not.direct) call runit(integral_package)
      call runit('xvmol2ja')

      if (mrcc.or.
     &    (iflags(18).ne.0.and.iflags(18).ne.3.or.iflags(23).ne.0.or.
     &     iflags(24).ne.0.or. iflags(25).ne.0.or.iflags(87).ne.0)
     &   ) call runit('xvprops')
      if (mrcc) call runit('rm -f VPOUT')
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
                  If (dkh .and. contract) Call Runit("xdkh")
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
         else
            if (.not.dirmp2) call runvcc
            if (iflags(87).ne.0) then
               call runit('xlambda')
               call runit('xvee')
            end if
         end if
      end if

      return
      end

