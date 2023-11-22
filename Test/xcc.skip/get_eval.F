
c This routine loads the SCF eigenvalues.

c INPUT
c    double precision dcore(*) : address space for xcc.com

c#define _PRINT_EIGENVALUES

      subroutine get_eval(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c INTERNAL VARIABLES
      integer irrep, isize
#ifdef _PRINT_EIGENVALUES
      integer i, ispin
#endif /* _PRINT_EIGENVALUES */

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"
#include "machsp.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c   o ALPHA spin
      do irrep = 1, nirrep
         ndx_eval(irrep,1,1) = freecore
         freecore = freecore + pop(irrep,1)
      end do
      do irrep = 1, nirrep
         ndx_eval(irrep,2,1) = freecore
         freecore = freecore + vrt(irrep,1)
      end do
      if (freecore.gt.1+dcoresize) then
         print *, '@GET_EVAL: Not enough memory.'
         call xcc_errex
      end if
      isize = freecore - ndx_eval(1,1,1)
      call getrec(20,'JOBARC','SCFEVALA',isize*iintfp,
     &            dcore(ndx_eval(1,1,1)))

c   o BETA spin
      if (iuhf.eq.0) then
         do irrep = 1, nirrep
            ndx_eval(irrep,1,2) = ndx_eval(irrep,1,1)
         end do
         do irrep = 1, nirrep
            ndx_eval(irrep,2,2) = ndx_eval(irrep,2,1)
         end do
      else
         do irrep = 1, nirrep
            ndx_eval(irrep,1,2) = freecore
            freecore = freecore + pop(irrep,2)
         end do
         do irrep = 1, nirrep
            ndx_eval(irrep,2,2) = freecore
            freecore = freecore + vrt(irrep,2)
         end do
         if (freecore.gt.1+dcoresize) then
            print *, '@GET_EVAL: Not enough memory.'
            call xcc_errex
         end if
         isize = freecore - ndx_eval(1,1,2)
         call getrec(20,'JOBARC','SCFEVALB',isize*iintfp,
     &               dcore(ndx_eval(1,1,2)))
      end if

#ifdef _PRINT_EIGENVALUES
      do ispin = 1, 1+iuhf
         do irrep = 1, nirrep
            isize = pop(irrep,ispin)
            print *, 'There are ',isize,' occupied orbitals in ',
     &               'irrep ',irrep,' with spin ',ispin,'.'
            do i = 0, isize-1
               print '(i5,a,i1,a,f16.10)',
     &            1+i,'[',irrep,'] : ',dcore(ndx_eval(irrep,1,ispin)+i)
            end do
         end do
         do irrep = 1, nirrep
            isize = vrt(irrep,ispin)
            print *, 'There are ',isize,' virtual orbitals in ',
     &               'irrep ',irrep,' with spin ',ispin,'.'
            do i = 0, isize-1
               print '(i5,a,i1,a,f16.10)',
     &            1+i,'[',irrep,'] : ',dcore(ndx_eval(irrep,2,ispin)+i)
            end do
         end do
      end do
#endif /* _PRINT_EIGENVALUES */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

