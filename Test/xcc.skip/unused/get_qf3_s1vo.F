
c This routine allocates space for and retrieves the 1-particle intermediate
c needed to compute the Qf energy from triples.

c INPUT
c    double precision dcore(*) : address space for xcc.com

#include "xcc.h"

c#define _PRINT_INTERMEDIATES

      subroutine get_qf3_s1vo(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c INTERNAL VARIABLES
      integer irrep, dim_scr, isize

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c   o ALPHA spin
      do irrep = 1, nirrep
         ndx_s1vo(irrep,1) = freecore
         freecore = freecore + (pop(irrep,1)*vrt(irrep,1))
      end do
      dim_scr = 1 + dcoresize - freecore
      if (dim_scr.lt.1) then
         write(*,*) '@GET_QF3_S1VO: Not enough memory.'
         call xcc_errex
      end if
      call make_qf3_s1vo(dcore,dcore(freecore),dim_scr,1)

c   o BETA spin
      if (iuhf.eq.0) then
         do irrep = 1, nirrep
            ndx_s1vo(irrep,2) = ndx_s1vo(irrep,1)
         end do
      else
         do irrep = 1, nirrep
            ndx_s1vo(irrep,2) = freecore
            freecore = freecore + (pop(irrep,2)*vrt(irrep,2))
         end do
         dim_scr = 1 + dcoresize - freecore
         if (dim_scr.lt.1) then
            write(*,*) '@GET_QF3_S1VO: Not enough memory.'
            call xcc_errex
         end if
         call make_qf3_s1vo(dcore,dcore(freecore),dim_scr,2)
      end if

#ifdef _PRINT_INTERMEDIATES
      do irrep = 1, nirrep
         isize = pop(irrep,1)*vrt(irrep,1)
         write(*,*) 'There are ',isize,' S1A amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_s1vo(irrep,1)),
     &                 vrt(irrep,1),pop(irrep,1),
     &                 vrt(irrep,1))
      end do
      if (iuhf.eq.1) then
      do irrep = 1, nirrep
         isize = pop(irrep,2)*vrt(irrep,2)
         write(*,*) 'There are ',isize,' S1B amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_s1vo(irrep,2)),
     &                 vrt(irrep,2),pop(irrep,2),
     &                 vrt(irrep,2))
      end do
      end if
#endif /* _PRINT_INTERMEDIATES */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

