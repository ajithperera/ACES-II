
c This routine tries to measure how much memory is required.

c OUTPUT
c    integer icoredim : the suggested size of icore to allocate
c    integer dcoredim : the number doubles after the silly I/O cache

#include "disttype.h"
#include "xcc.h"

#define BASE_DCOREDIM 1024

      subroutine xcc_measure(icoredim,dcoredim)
      implicit none

c ARGUMENT LIST
      integer icoredim, dcoredim

c INTERNAL VARIABLES
c   o maximize various dimensions
      integer max_pop, max_vrt, max_irpdpd(22)
      integer max_oo, max_vv, max_ov
c   o maximize the amount of scratch space
      integer max_scr, max_tmp, t3_scr
c   o count T3 amplitudes
      integer irp_abc, t3occ(8), t3vrt(8)
c   o other junk
      integer file_cache
      integer i, ispin, itmp, irrep

c COMMON BLOCKS
#include "flags.com"
#include "machsp.com"
#include "method.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "info.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _DEBUG_XCC
      print *, '\n @XCC: ---=== MEMORY MEASUREMENTS ===---'
#endif

      max_scr = 0

c   o maximize various dimensions

      max_pop = 0
      do ispin = 1, 1+iuhf
         do irrep = 1, nirrep
            max_pop = max(max_pop,pop(irrep,ispin))
         end do
      end do

      max_vrt = 0
      do ispin = 1, 1+iuhf
         do irrep = 1, nirrep
            max_vrt = max(max_vrt,vrt(irrep,ispin))
         end do
      end do

      do i = 1, 22
         max_irpdpd(i) = 0
         do irrep = 1, nirrep
            max_irpdpd(i) = max(max_irpdpd(i),irpdpd(irrep,i))
         end do
      end do
      max_oo = max(max_irpdpd(h_DistType_IltJ),
     &             max_irpdpd(h_DistType_iltj),
     &             max_irpdpd(h_DistType_Ij)   )
      max_vv = max(max_irpdpd(h_DistType_AltB),
     &             max_irpdpd(h_DistType_altb),
     &             max_irpdpd(h_DistType_Ab)   )
      max_ov = max(max_irpdpd(h_DistType_AI),
     &             max_irpdpd(h_DistType_ai),
     &             max_irpdpd(h_DistType_Ai),
     &             max_irpdpd(h_DistType_Ia) )

c ----------------------------------------------------------------------

c   o What's a few floats among friends, eh?
      dcoredim = BASE_DCOREDIM

c   o SCF Eigenvalues
      itmp = nocco(1) + nvrto(1)
      if (iuhf.eq.1) itmp = itmp + nocco(2) + nvrto(2)
      dcoredim = dcoredim + itmp
#ifdef _DEBUG_XCC
      print *, 'SCF eigenvalues require: ',(itmp*ifltln),' Bytes'
#endif

c   o T1
      itmp = irpdpd(1,h_DistType_AI)
      if (iuhf.eq.1) itmp = itmp + irpdpd(1,h_DistType_ai)
      dcoredim = dcoredim + itmp
#ifdef _DEBUG_XCC
      print *, 'T1 amplitudes require: ',(itmp*ifltln),' Bytes'
#endif

c   o T2 and T2(1)
      itmp = 0
      do irrep = 1, nirrep
         itmp =   itmp
     &          + (   irpdpd(irrep,h_DistType_AltB)
     &              * irpdpd(irrep,h_DistType_IltJ) )
      end do
      if (iuhf.eq.1) then
         do irrep = 1, nirrep
            itmp =   itmp
     &             + (   irpdpd(irrep,h_DistType_altb)
     &                 * irpdpd(irrep,h_DistType_iltj) )
         end do
      end if
      do irrep = 1, nirrep
         itmp =   itmp
     &          + (   irpdpd(irrep,h_DistType_Ab)
     &              * irpdpd(irrep,h_DistType_Ij) )
      end do
      dcoredim = dcoredim + ( 2 * itmp )
#ifdef _DEBUG_XCC
      print *, 'T2 amplitudes require: ',(itmp*ifltln),' Bytes'
      print *, 'T2(1) amplitudes require: ',(itmp*ifltln),' Bytes'
#endif

c   o adj(T2)*T1^2 scratch vector
      itmp = max_irpdpd(h_DistType_AltB)
      if (iuhf.eq.1) itmp = max(itmp,max_irpdpd(h_DistType_altb))
      itmp = max(itmp,max_vrt)
      max_scr = max(max_scr,itmp)

c   o adj(T1)^2*T1^2 scratch vector
      itmp = 2*max_vrt*max_vrt
      max_scr = max(max_scr,itmp)

c ----------------------------------------------------------------------

      t3_scr = 0

#ifdef ONLY_ONE_T3_DIST

c   o T3AAB & T3BBA
      itmp = 0
      do ispin = 1, 1+iuhf
         call irpdpd_pltpq(vrt(1,ispin),vrt(1,3-ispin),t3vrt)
         do irp_abc = 1, nirrep
            itmp = max(itmp,t3vrt(irp_abc))
         end do
c     end do ispin = 1, 1+iuhf
      end do
      dcoredim = dcoredim + itmp
#ifdef _DEBUG_XCC
      print *, 'The largest T3 distribution requires: ',itmp*ifltln,
     &         ' Bytes'
#endif /* _DEBUG_XCC */

#else /* store all the T3 amplitudes */

c   o T3AAA & T3BBB
      do ispin = 1, 1+iuhf
         call irpdpd_pltpltp(pop(1,ispin),t3occ)
         call irpdpd_pltpltp(vrt(1,ispin),t3vrt)
         itmp = 0
         do irp_abc = 1, nirrep
            itmp =   itmp
     &             + ( t3occ(irp_abc) * t3vrt(irp_abc) )
         end do
         dcoredim = dcoredim + itmp
#ifdef _DEBUG_XCC
         if (ispin.eq.1) then
            print *, 'T3AAA amplitudes require: ',itmp*ifltln,' Bytes'
         else
            print *, 'T3BBB amplitudes require: ',itmp*ifltln,' Bytes'
         end if
#endif /* _DEBUG_XCC */
c     end do ispin = 1, 1+iuhf
      end do

c   o T3AAB & T3BBA
      do ispin = 1, 1+iuhf
         call irpdpd_pltpq(pop(1,ispin),pop(1,3-ispin),t3occ)
         call irpdpd_pltpq(vrt(1,ispin),vrt(1,3-ispin),t3vrt)
         itmp = 0
         do irp_abc = 1, nirrep
            itmp =   itmp
     &             + ( t3occ(irp_abc) * t3vrt(irp_abc) )
         end do
         dcoredim = dcoredim + itmp
#ifdef _DEBUG_XCC
         if (ispin.eq.1) then
            print *, 'T3AAB amplitudes require: ',itmp*ifltln,' Bytes'
         else
            print *, 'T3BBA amplitudes require: ',itmp*ifltln,' Bytes'
         end if
#endif /* _DEBUG_XCC */
c     end do ispin = 1, 1+iuhf
      end do

#endif /* ONLY_ONE_T3_DIST */

c   o scratch space for one alpha/beta/mixed T3(ijk,pqr) product vector
      itmp = 0
      do ispin = 1, 1+iuhf
         call irpdpd_pltpltp(vrt(1,ispin),t3vrt)
         do irp_abc = 1, nirrep
            max_tmp = max(t3vrt(irp_abc),itmp)
         end do
      end do
      do ispin = 1, 1+iuhf
         call irpdpd_pltpq(vrt(1,ispin),vrt(1,3-ispin),t3vrt)
         do irp_abc = 1, nirrep
            max_tmp = max(t3vrt(irp_abc),itmp)
         end do
      end do
      t3_scr = t3_scr + max_tmp

c   o scratch space for T2(mi,pq) or T2(jk,er) dgemm factors
      itmp = 0
      do i = 1, nirrep
         itmp = itmp + ( pop(i,1) * max_vv )
      end do
      max_tmp = itmp
      itmp = 0
      do i = 1, nirrep
         itmp = itmp + ( pop(i,2) * max_vv )
      end do
      max_tmp = max(max_tmp,itmp)
      max_tmp = max(max_tmp,max_vrt)
      t3_scr = t3_scr + max_tmp

c   o scratch space for one alpha/beta/mixed 2-electron (ab,ci) integral
      t3_scr = t3_scr + max_vv

c   o scratch space for all alpha/beta/mixed 2-electron (ij,ka) integrals
c     For calculating T3AAB/BBA, we need all three sets:
c        <Ij|Ka>, <iJ|kA>, and <IJ||KA>
c     in core this is because I am lazy and I don't want to reload these
c     lists over and over for each I, J, k combo.
      itmp = 0
      do i = 1, nirrep
         itmp =   itmp
     &          + (   irpdpd(i,h_DistType_IltJ)
     &              * irpdpd(i,h_DistType_AI)   )
      end do
      do i = 1, nirrep
         itmp =   itmp
     &          + (   irpdpd(i,h_DistType_Ij)
     &              * irpdpd(i,h_DistType_Ia) )
      end do
      if (iuhf.eq.1) then
         do i = 1, nirrep
            itmp =   itmp
     &             + (   irpdpd(i,h_DistType_Ij)
     &                 * irpdpd(i,h_DistType_Ai) )
         end do
      end if
      max_tmp = itmp
      itmp = 0
      do i = 1, nirrep
         itmp =   itmp
     &          + (   irpdpd(i,h_DistType_iltj)
     &              * irpdpd(i,h_DistType_ai)   )
      end do
      do i = 1, nirrep
         itmp =   itmp
     &          + (   irpdpd(i,h_DistType_Ij)
     &              * irpdpd(i,h_DistType_Ai) )
      end do
      if (iuhf.eq.1) then
         do i = 1, nirrep
            itmp =   itmp
     &             + (   irpdpd(i,h_DistType_Ij)
     &                 * irpdpd(i,h_DistType_Ia) )
         end do
      end if
      max_tmp = max(max_tmp,itmp)
      t3_scr = t3_scr + max_tmp

c   o add the T3 scratch space to the total scratch space
      max_scr = max(max_scr,t3_scr)
#ifdef _DEBUG_XCC
      print *, 'T3 scratch needs: ',t3_scr,' floats ( = ',
     &         t3_scr*ifltln,' Bytes )'
#endif

c ----------------------------------------------------------------------

c   o collect the scratch space
      dcoredim = dcoredim + max_scr
#ifdef _DEBUG_XCC
      print *, 'total scratch needs: ',(max_scr*ifltln),' Bytes'
#endif

c   o Finish with the silly file cache. (cf. libr/popcorn.F)
      file_cache = iflags(37)*iflags(44)
#ifdef _DEBUG_XCC
      print *, 'file cache: ',(file_cache*iintln),' Bytes'
#endif

c   o Get it while it's good...
      itmp = -1 + ( ( iflags(36) - file_cache ) / iintfp )
      dcoredim = max(itmp,dcoredim)

c   o Close the deal.
      icoredim = file_cache + ( dcoredim * iintfp )
      icoredim = icoredim + mod(icoredim,iintfp)
#ifdef _DEBUG_XCC
      itmp = icoredim*iintln
      itmp = rshift(itmp,20) + 1
      print *, '\n          TOTAL MEMORY NEEDS: ',itmp,' MB'
      print *, '      ---===========================---\n'
#endif

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

