      subroutine symcntr(atom,coord,integaxis,numoct)

c This routine determines the integration requirement for each of the
c axis around the current atom based on symmetry.

      implicit none
#include <aces.h>
#include <mol.com>
      integer atom,integaxis(3,natoms),numoct
      M_REAL coord(3,natoms)

      integer ixyz,nzero,first0,first

      call callstack_push('SYMCNTR')

      numoct=8

c nzero  : the number (0-3) of coordinates which are equal to zero
c first0 : the first zero coordinate
c first  : the first non-zero coordinate

      nzero=0
      first0=4
      first=4
      do ixyz=3,1,-1
        integaxis(ixyz,atom)=0
        if (coord(ixyz,atom).eq.0.d0) nzero=nzero+1
        if (coord(ixyz,atom).eq.0.d0) first0=ixyz
        if (coord(ixyz,atom).ne.0.d0) first=ixyz
      end do
c Determine which axis need to be fully integrated over and which can
c be integrated over the positive part only.  This is based on which
c cartesian coordinates of the atom are zero.
c
c numoct    The number of types of octants.
c integaxis If set for a given axis, an integral needs to only integrate
c           over positive values (i.e. for the D2h case with x and y zero,
c           integaxi=1,2,0 - this means that you have to integrate over
c           all z but only positive x and y)

      if (compptgp.eq.'D2h ') then

c       zero   numoct  axis
c
c       x y z    1    1  2  3
c
c       x y -    2    1  2  -
c       x - z    2    1  -  3
c       - y z    2    -  2  3
c
c       x - -    4    1  -  -
c       - y -    4    -  2  -
c       - - z    4    -  -  3
c
c       - - -    8    -  -  -

        if      (nzero.eq.3) then
          numoct=1
          integaxis(1,atom)=1
          integaxis(2,atom)=2
          integaxis(3,atom)=3
        else if (nzero.eq.2) then
          numoct=2
          integaxis(1,atom)=1
          integaxis(2,atom)=2
          integaxis(3,atom)=3
          integaxis(first,atom)=0
        else if (nzero.eq.1) then
          numoct=4
          integaxis(first0,atom)=first0
        endif

      else if (compptgp.eq.'C2h ') then

c       zero   numoct  axis
c
c       x y z    2    1  -  3  (either one)
c                     -  2  3
c
c       ? ? z    4    -  -  3
c       x y -    4    1  -  -  (either one)
c                     -  2  -
c
c       - - -    8    -  -  -

        if      (nzero.eq.3) then
          numoct=2
          integaxis(1,atom)=1
          integaxis(3,atom)=3
        else if (coord(3,atom).eq.0.d0) then
          numoct=4
          integaxis(3,atom)=3
        else if (nzero.eq.2) then
          numoct=4
          integaxis(1,atom)=1
        endif

      else if (compptgp.eq.'D2  ') then

c       Fixed by SI 12/10/98. There are
c       more problems in cmp2zmat transformation.

c       zero   numoct  axis
c
c       x y z    2    1  2  -  (either one)
c                     1  -  3
c                     -  2  3
c
c       x y -    4    -  -  3
c       x - z    4    -  2  -
c       - y z    4    1  -  -
c
c       ? ? ?    8    -  -  -

        if      (nzero.eq.3) then
          numoct=2
          integaxis(1,atom)=1
          integaxis(2,atom)=2
        else if (nzero.eq.2) then
          numoct=4
          integaxis(first0,atom)=first0
        endif

      else if (compptgp.eq.'C2v ') then

c       zero   numoct  axis
c
c       x y ?    2    1  2  -
c
c       x - ?    4    1  -  -
c       - y ?    4    -  2  -
c
c       - - ?    8    -  -  -

        if (first.gt.2) then
          numoct=2
          integaxis(1,atom)=1
          integaxis(2,atom)=2
        else if (first0.le.2) then
          numoct=4
          integaxis(first0,atom)=first0
        endif

      else if (compptgp.eq.'C2  ') then

c       zero   numoct  axis
c
c       x y ?    4    1  -  -  (either one)
c                     -  2  -
c
c       ? ? ?    8    -  -  -

        if (first.gt.2) then
          numoct=4
          integaxis(1,atom)=1
        endif

      else if (compptgp.eq.'C i ') then

c       zero   numoct  axis
c
c       x y z    4    1  -  -  (either one)
c                     -  2  -
c                     -  -  3
c
c       ? ? ?    8    -  -  -

        if (nzero.eq.3) then
          numoct=4
          integaxis(1,atom)=1
        endif

      else if (compptgp.eq.'C s ')then

c       zero   numoct  axis
c
c       ? ? z    4    -  -  3
c
c       ? ? ?    8    -  -  -

        if (coord(3,atom).eq.0.d0) then
          numoct=4
          integaxis(3,atom)=3
        endif

      endif

      call callstack_pop
      return
      end
