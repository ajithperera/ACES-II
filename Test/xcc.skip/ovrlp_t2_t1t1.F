
c This routine adds the quantity <0|adj(T2)*T1^2|0>/2 and its complex
c conjugate to the overlap.

c INPUT
c    double precision dcore(*) : address space for xcc.com

      subroutine ovrlp_t2_t1t1(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c INTERNAL VARIABLES
      double precision ovrlp_tmp

c COMMON BLOCKS
#include "xcc.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      ovrlp_tmp = 0.0d0
      call ovrlp_t2aa_t1t1(dcore,dcore(freecore),ovrlp_tmp)
      call ovrlp_t2ab_t1t1(dcore,dcore(freecore),ovrlp_tmp)

      overlap = overlap + ovrlp_tmp
      print *, '<0|adj(T2)*T1^2|0>/2 + c.c. = ',ovrlp_tmp

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

