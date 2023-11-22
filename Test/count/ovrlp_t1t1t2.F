
c This routine is a driver for the subroutines that add the quantity
c <0|adj(T2)*T1^2|0> and its complex conjugate to the overlap. It
c mainly maps the address of dcore(freecore) into dscratch(1).

c INPUT
c    double precision dcore(*) : core space for doubles

      subroutine ovrlp_t1t1t2(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "syminf.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c    o add the T2AA/T2BB piece
       if (nirrep.eq.1) then
          call nosym(dcore,1)
          if (iuhf.eq.1) call nosym(dcore,2)
       else
          call ovrlp_t1t1t2aa(dcore,dcore(freecore),1)
          if (iuhf.eq.1) call ovrlp_t1t1t2aa(dcore,dcore(freecore),2)
       end if

c    o add the T2AB piece
c       call ovrlp_t1t1t2ab(dcore,dcore(freecore))

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

