
c This routine loads the T1 amplitudes and calculates <0|adj(T1)*T1|0>.

c INPUT
c    double precision dcore(*) : address space for xcc.com

#include "disttype.h"
#include "listnum.h"

c#define _PRINT_AMPLITUDES

      subroutine get_t1(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c EXTERNAL FUNCTIONS
      double precision ddot, dnormsqr

c INTERNAL VARIABLES
      integer irrep, isize
      integer listnum, sublist
      double precision dtmp

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      t_sqr(1) = 0.0d0

c   o ALPHA spin
      listnum = h_T1_IA_AI
      sublist = 1
      do irrep = 1, nirrep
         ndx_t1(irrep,1) = freecore
         freecore = freecore + (pop(irrep,1)*vrt(irrep,1))
      end do
      if (freecore.gt.1+dcoresize) then
         print *, '@GET_T1: Not enough memory.'
         call xcc_errex
      end if
      isize = irpdpd(1,h_DistType_AI)
      call getlst(dcore(ndx_t1(1,1)),1,1,1,sublist,listnum)
c      dtmp = ddot(isize,dcore(ndx_t1(1,1)),1,
c     &                  dcore(ndx_t1(1,1)),1)
      dtmp = dnormsqr(isize,dcore(ndx_t1(1,1)),1)
      t_sqr(1) = dtmp
#ifdef _DEBUG_XCC
      print *, '<0|adj(T1A)*T1A|0> = ', dtmp
#endif /* _DEBUG_XCC */

c   o BETA spin
      if (iuhf.eq.0) then
         t_sqr(1) = t_sqr(1) + dtmp
         do irrep = 1, nirrep
            ndx_t1(irrep,2) = ndx_t1(irrep,1)
         end do
      else
         listnum = h_T1_ia_ai
         sublist = 2
         do irrep = 1, nirrep
            ndx_t1(irrep,2) = freecore
            freecore = freecore + (pop(irrep,2)*vrt(irrep,2))
         end do
         if (freecore.gt.1+dcoresize) then
            print *, '@GET_T1: Not enough memory.'
            call xcc_errex
         end if
         isize = irpdpd(1,h_DistType_ai)
         call getlst(dcore(ndx_t1(1,2)),1,1,1,sublist,listnum)
c         dtmp = ddot(isize,dcore(ndx_t1(1,2)),1,
c     &                     dcore(ndx_t1(1,2)),1)
         dtmp = dnormsqr(isize,dcore(ndx_t1(1,2)),1)
         t_sqr(1) = t_sqr(1) + dtmp
      end if
#ifdef _DEBUG_XCC
      print *, '<0|adj(T1B)*T1B|0> = ', dtmp
#endif

      overlap = overlap + t_sqr(1)
      print *, '<0|adj(T1)*T1|0> = ', t_sqr(1)

#ifdef _PRINT_AMPLITUDES
      do irrep = 1, nirrep
         isize = pop(irrep,1)*vrt(irrep,1)
         print *, 'There are ',isize,' T1A amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_t1(irrep,1)),
     &                 vrt(irrep,1),pop(irrep,1),
     &                 vrt(irrep,1))
      end do
      if (iuhf.eq.1) then
      do irrep = 1, nirrep
         isize = pop(irrep,2)*vrt(irrep,2)
         print *, 'There are ',isize,' T1B amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_t1(irrep,2)),
     &                 vrt(irrep,2),pop(irrep,2),
     &                 vrt(irrep,2))
      end do
      end if
#endif /* _PRINT_AMPLITUDES */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

