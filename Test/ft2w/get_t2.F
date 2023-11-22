
c This routine loads the T2 amplitudes.

c INPUT
c double dcore(*) : address space for ft2w.com

#include "disttype.h"
#include "listnum.h"

c#define _PRINT_AMPLITUDES

      subroutine get_t2(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c INTERNAL VARIABLES
      integer iRight, i, j
      integer irrep, ioff, itmp

c COMMON BLOCKS
#include "ft2w.com"
#include "method.com" /* for iuhf */
#include "syminf.com" /* for nirrep */
#include "sympop.com" /* for irpdpd(,) */

c ----------------------------------------------------------------------

c   o T2(IJ,AB); A<B, I<J
      iRight = h_T2_IJAB_AltB_IltJ
      do irrep = 1, nirrep
         ndx_t2(irrep,1) = freecore
         i = irpdpd(irrep,h_DistType_AltB)
         j = irpdpd(irrep,h_DistType_IltJ)
         ioff = i * j
         if (ioff.ne.0) then
            itmp = freecore + ioff
            if (itmp.gt.dcoresize) then
               print *, '@GET_T2: Not enough memory.'
               call aces_exit(1)
            end if
            call getlst(dcore(ndx_t2(irrep,1)),1,j,2,irrep,iRight)
            freecore = itmp
         end if
      end do

c   o T2(ij,ab); a<b, i<j
      if (iuhf.eq.0) then
         do irrep = 1, nirrep
            ndx_t2(irrep,2) = ndx_t2(irrep,1)
         end do
      else
         iRight = h_T2_ijab_altb_iltj
         do irrep = 1, nirrep
            ndx_t2(irrep,2) = freecore
            i = irpdpd(irrep,h_DistType_altb)
            j = irpdpd(irrep,h_DistType_iltj)
            ioff = i * j
            if (ioff.ne.0) then
               itmp = freecore + ioff
               if (itmp.gt.dcoresize) then
                  print *, '@GET_T2: Not enough memory.'
                  call aces_exit(1)
               end if
               call getlst(dcore(ndx_t2(irrep,2)),1,j,2,irrep,iRight)
               freecore = itmp
            end if
         end do
      end if

c   o T2(Ij,Ab); Ab, Ij
      iRight = h_T2_IjAb_Ab_Ij
      do irrep = 1, nirrep
         ndx_t2(irrep,3) = freecore
         i = irpdpd(irrep,h_DistType_Ab)
         j = irpdpd(irrep,h_DistType_Ij)
         ioff = i * j
         if (ioff.ne.0) then
            itmp = freecore + ioff
            if (itmp.gt.dcoresize) then
               print *, '@GET_T2: Not enough memory.'
               call aces_exit(1)
            end if
            call getlst(dcore(ndx_t2(irrep,3)),1,j,2,irrep,iRight)
            freecore = itmp
         end if
      end do

c ----------------------------------------------------------------------

#ifdef _PRINT_AMPLITUDES
      do irrep = 1, nirrep
         i = irpdpd(irrep,h_DistType_AltB)
         j = irpdpd(irrep,h_DistType_IltJ)
         print *, 'There are ',i*j,' T2AA amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_t2(irrep,1)),i,j,i)
      end do
      if (iuhf.eq.1) then
      do irrep = 1, nirrep
         i = irpdpd(irrep,h_DistType_altb)
         j = irpdpd(irrep,h_DistType_iltj)
         print *, 'There are ',i*j,' T2BB amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_t2(irrep,2)),i,j,i)
      end do
      end if
      do irrep = 1, nirrep
         i = irpdpd(irrep,h_DistType_Ab)
         j = irpdpd(irrep,h_DistType_Ij)
         print *, 'There are ',i*j,' T2AB amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_t2(irrep,3)),i,j,i)
      end do
#endif /* _PRINT_AMPLITUDES */

c ----------------------------------------------------------------------

      return
c     end subroutine get_t2
      end

