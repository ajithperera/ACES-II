
c This routine tries to measure how much memory is required.

c OUTPUT
c    integer icoredim : the suggested size of icore to allocate
c    integer dcoredim : the number doubles after the silly I/O cache

#include "disttype.h"

      subroutine xcc_measure(icoredim,dcoredim)
      implicit none

c ARGUMENT LIST
      integer icoredim, dcoredim

c INTERNAL VARIABLES
c   o maximize the amount of scratch space
      integer scr_max
      integer i, itmp

c COMMON BLOCKS
#include "flags.com"
#include "machsp.com"
#include "method.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c   o What's a few floats among friends, eh?
      dcoredim=1024

c   o T1
      itmp=irpdpd(1,h_DistType_AI)
      if (iuhf.eq.1) then
         itmp=itmp+irpdpd(1,h_DistType_ai)
      end if
      dcoredim=dcoredim+itmp

c   o T2
      itmp=0
      do i = 1,nirrep
         itmp=itmp+irpdpd(i,h_DistType_AltB)*irpdpd(i,h_DistType_IltJ)
      end do
      if (iuhf.eq.1) then
         do i = 1,nirrep
            itmp=itmp+irpdpd(i,h_DistType_altb)*irpdpd(i,h_DistType_iltj)
         end do
      end if
      do i = 1,nirrep
         itmp=itmp+irpdpd(i,h_DistType_Ab)*irpdpd(i,h_DistType_Ij)
      end do
      dcoredim=dcoredim+itmp

c   o adj(T2)*T1^2 scratch vector
      itmp=0
      do i = 1,nirrep
         itmp=max(itmp,irpdpd(i,h_DistType_AltB))
      end do
      if (iuhf.eq.1) then
         do i = 1,nirrep
            itmp=max(itmp,irpdpd(i,h_DistType_altb))
         end do
      end if
      do i = 1,nirrep
         itmp=max(itmp,vrt(i,1))
      end do
      scr_max=max(itmp,scr_max)

c   o collect the scratch space
      dcoredim = dcoredim + scr_max

c   o Finish with the silly file cache. (cf. libr/popcorn.F)
      itmp     = iflags(37)*iflags(44)
      icoredim = itmp + ( dcoredim * iintfp )
      icoredim = icoredim + mod(icoredim,iintfp)

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

