
c This routine is the C1 equivalent of ovrlp_t1t1t2aa.F. It adds the
c quantity <0|adj(T2AA)*T1A^2|0> and its complex conjugate to the overlap.
c It can also be used for the T2BB piece if the argument ispin is 2;
c otherwise, if iuhf is 0, it doubles the T2AA term.

c INPUT
c    double precision dcore(*) : core space for doubles
c    integer          ispin    : Alpha/Beta signal (1=Alpha,2=Beta)

#include "disttype.h"

      subroutine nosym(dcore,ispin)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)
      integer ispin

c INTERNAL VARIABLES
      double precision dtmp, ovrlp_tmp
      integer t2_start, t1_start, ioff, nocc, nvrt
      integer i, j, a, b

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      nocc = pop(1,ispin)
      if (nocc.lt.2) return
      nvrt = vrt(1,ispin)
      if (nvrt.lt.2) return

      ioff = 0
      ovrlp_tmp = 0.0d0
      t1_start = ndx_t1(1,ispin)
      t2_start = ndx_t2(1,ispin)

      do j = 1,nocc-1
      do i = 0,j-1
      do b = 1,nvrt-1
         dtmp = dcore(t1_start+nvrt*j+b)
      do a = 0,b-1
         write(*,'(3F16.10,3x,5I5)')
     &      dcore(t2_start+ioff),
     &      dtmp,
     &      dcore(t1_start+nvrt*i+a),
     &      t2_start+ioff, t1_start+nvrt*j+b, t1_start+nvrt*i+a, a, i
c     &      a+1, b+1, i+1, j+1
         ovrlp_tmp =   ovrlp_tmp
     &               + (   dtmp
     &                   * dcore(t1_start+nvrt*i+a)
     &                   * dcore(t2_start+ioff) )
         ioff = ioff + 1
      end do
      end do
      end do
      end do

      ovrlp_tmp = 2.0d0*ovrlp_tmp
      if (iuhf.eq.0) then
         overlap = overlap + ovrlp_tmp*2.0d0
         write(*,*) '@XCC: <0|adj(T2AA)*T1A*T1A|0> = ',ovrlp_tmp
         write(*,*) '@XCC: <0|adj(T2BB)*T1B*T1B|0> = ',ovrlp_tmp
      else
         overlap = overlap + ovrlp_tmp
         if (ispin.eq.1) then
            write(*,*) '@XCC: <0|adj(T2AA)*T1A*T1A|0> = ',ovrlp_tmp
         else
            write(*,*) '@XCC: <0|adj(T2BB)*T1B*T1B|0> = ',ovrlp_tmp
         end if
      end if

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

