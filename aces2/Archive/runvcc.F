
      subroutine runvcc

#include "aces2.h"

#ifdef _TRACE_ACES2
      write(*,*) '@ACES2: Entered ', __FILE__
#endif
      character*79 szEExtrap
      szEExtrap = 'xa2proc ener_extrp'
c
c KJW 9/29/99
c support for Stan Kucharski's' ccsdtq program
c
      if (iflags2(152).eq.1) then
         call runit('xrdycom')
         call runit('xccsdtq < COM')
         return
      end if

1     call runit('xvcc')

c JDW 3/31/94
c noniterative 5th-order triples and quadruples
c
    if (iflags(2).eq.12.or.(iflags(2).ge.26.and.iflags(2).le.30)) then
         call runit('xvcc5t')
         call runit('xvcc5q')
      end if

c JDW 3/17/95
c CC5SD[T]
      if (iflags(2).eq.31) call runit('xvcc5t')
c
c Extrapolate energy (only for correlated methods), 01/2006, Ajith Perera

      if (bEExtrap .or. bCExtrap) then
          call runit(szEExtrap)
      endif
c
      if (iflags(22).eq.1) then
         call a2getrec(20,'JOBARC','BRUKTEST',1,itest)
         if(itest.eq.1)return
         if (nddo_guess) call runit('xnddo')
         call runit('xvscf')
         call runit('xvtran')
         call runit('xintprc')
         goto 1
      end if

      end

