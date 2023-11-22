
      subroutine runoombpt
      implicit none
#include "aces2.h"
      integer ndrgeo,itest
      logical bDoneOO

      bDoneOO = .false.
      ndrgeo = 0
      itest = 0

      call a2getrec(20,'JOBARC','IFLAGS',100,iflags)
      do while (.not. bDoneOO)

        call runit('xoombpt2')
        call a2getrec(20,'JOBARC','IFLAGS',100,iflags)
        call a2getrec(20,'JOBARC','OOMBPTST',1,itest)

C If itest == 2 then error
        if (itest .eq. 2) stop
C If itest == 0 or 1 run xvscf and precc
        call a2putrec(20,'JOBARC','NOTRANS ',1,0)
        If (dkh .and. contract) Call Runit("xdkh")
        call runit('xvscf')
        call precc(ndrgeo)

        if (itest .eq. 1) bDoneOO = .true.

      end do

      iflags(2) = 1
      call a2putrec(20,'JOBARC','IFLAGS',100,iflags)

      call runvcc

      return
      end
