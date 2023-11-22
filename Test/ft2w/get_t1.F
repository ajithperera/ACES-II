
c This routine loads the T1 amplitudes.

c INPUT
c double dcore(*) : address space for ft2w.com

#include "disttype.h"
#include "listnum.h"

c#define _PRINT_AMPLITUDES

      subroutine get_t1(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c INTERNAL VARIABLES
      integer irrep, iSize
      integer iLeft, iRight

c COMMON BLOCKS
#include "ft2w.com"
#include "method.com" /* for iuhf */
#include "sym.com" /* for pop(,) and vrt(,) */
#include "syminf.com" /* for nirrep */
#include "sympop.com" /* for irpdpd(,) */

c ----------------------------------------------------------------------

c   o ALPHA spin
      iRight = h_T1_IA_AI
      iLeft  = 1
      do irrep = 1, nirrep
         ndx_t1(irrep,1) = freecore
         freecore = freecore + (vrt(irrep,1)*pop(irrep,1))
      end do
      if (freecore.gt.dcoresize) then
         print *, '@GET_T1: Not enough memory.'
         call aces_exit(1)
      end if
      iSize = irpdpd(1,h_DistType_AI)
      call getlst(dcore(ndx_t1(1,1)),1,1,1,iLeft,iRight)

c   o BETA spin
      if (iuhf.eq.0) then
         do irrep = 1, nirrep
            ndx_t1(irrep,2) = ndx_t1(irrep,1)
         end do
      else
         iRight = h_T1_ia_ai
         iLeft  = 2
         do irrep = 1, nirrep
            ndx_t1(irrep,2) = freecore
            freecore = freecore + (vrt(irrep,2)*pop(irrep,2))
         end do
         if (freecore.gt.dcoresize) then
            print *, '@GET_T1: Not enough memory.'
            call aces_exit(1)
         end if
         iSize = irpdpd(1,h_DistType_ai)
         call getlst(dcore(ndx_t1(1,2)),1,1,1,iLeft,iRight)
      end if

c ----------------------------------------------------------------------

#ifdef _PRINT_AMPLITUDES
      do irrep = 1, nirrep
         iSize = vrt(irrep,1)*pop(irrep,1)
         print *, 'There are ',iSize,' T1A amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_t1(irrep,1)),
     &                 vrt(irrep,1),pop(irrep,1),
     &                 vrt(irrep,1))
      end do
      if (iuhf.eq.1) then
      do irrep = 1, nirrep
         iSize = vrt(irrep,2)*pop(irrep,2)
         print *, 'There are ',iSize,' T1B amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_t1(irrep,2)),
     &                 vrt(irrep,2),pop(irrep,2),
     &                 vrt(irrep,2))
      end do
      end if
#endif /* _PRINT_AMPLITUDES */

c ----------------------------------------------------------------------

      return
c     end subroutine get_t1
      end

