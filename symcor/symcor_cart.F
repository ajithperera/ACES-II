      subroutine symcor_cart(refq, Rint, scr,nsize,nmode)
c
c Generate cartesian coordinates from a set of internals
c   generated in setpts4_pes_normal in symcor.
c This subroutine should be in symcor, but we use common blocks 
c that are used in joda exclusively. 
c The internals in Rint are overwrtten by cartesians
c The cartesians in Rint are aligned in maximum coincidence with refq
c
#include "mxatms.par"
#include "coord.com"
C
      double precision Refq(nsize), Rint(nsize), scr(*)
      logical print
c
      nrx = nsize
      nx6m = nmode
      print = .false.
c
      if (print) then
      write(6,*) '@symcor_cart: Rint on input, natoms', natoms
      call output(Rint, 1, 1, 1, nsize, 1, nsize, 1)
      endif
c
      call usqush(Rint,nx6m)
      call DCOPY(nrx, Rint, 1, R, 1)
      call gen_cart_coord_b(scr, print)
      call mn_align_geom(natoms, atmass, q, refq, scr, mxcor)
      call dcopy(nrx, q, 1, Rint, 1)
c
      if (print) then
      write(6,*) '@symcor_cart: Rint on output'
      call output(Rint, 1, 1, 1, nsize, 1, nsize, 1)
      endif

c
      return
      end
      
