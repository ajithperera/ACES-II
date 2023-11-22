
c This routine loads the T2 amplitudes and calculates <0|adj(T2)*T2|0>.

c INPUT
c    double precision dcore(*) : address space for xcc.com

#include "disttype.h"
#include "listnum.h"

c#define _PRINT_AMPLITUDES

      subroutine get_t2(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c EXTERNAL FUNCTIONS
      double precision ddot, dnormsqr

c INTERNAL VARIABLES
      integer listnum, i, j
      integer irrep, ioff, itmp
      double precision dtmp

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "machsp.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      t_sqr(2) = 0.0d0

c   o T2(IJ,AB); A<B, I<J
      listnum = h_T2_IJAB_AltB_IltJ
      do irrep = 1, nirrep
         ndx_t2(irrep,1) = freecore
         i = irpdpd(irrep,h_DistType_AltB)
         j = irpdpd(irrep,h_DistType_IltJ)
         ioff = i * j
         if (ioff.ne.0) then
            itmp = freecore + ioff
            if (itmp.gt.(1+dcoresize)) then
               print *, '@GET_T2: Not enough memory.'
               call xcc_errex
            end if
            call getlst(dcore(ndx_t2(irrep,1)),1,j,2,irrep,listnum)
            freecore = itmp
         end if
      end do
      itmp = freecore - ndx_t2(1,1)
c      dtmp = ddot(itmp,dcore(ndx_t2(1,1)),1,
c     &                 dcore(ndx_t2(1,1)),1)
      dtmp = dnormsqr(itmp,dcore(ndx_t2(1,1)),1)
      t_sqr(2) = t_sqr(2) + dtmp
#ifdef _DEBUG_XCC
      print *, '<0|adj(T2AA)*T2AA|0> = ', dtmp
#endif /* _DEBUG_XCC */

c   o T2(ij,ab); a<b, i<j
      if (iuhf.eq.0) then
         t_sqr(2) = t_sqr(2) + dtmp
         do irrep = 1, nirrep
            ndx_t2(irrep,2) = ndx_t2(irrep,1)
         end do
      else
         listnum = h_T2_ijab_altb_iltj
         do irrep = 1, nirrep
            ndx_t2(irrep,2) = freecore
            i = irpdpd(irrep,h_DistType_altb)
            j = irpdpd(irrep,h_DistType_iltj)
            ioff = i * j
            if (ioff.ne.0) then
               itmp = freecore + ioff
               if (itmp.gt.(1+dcoresize)) then
                  print *, '@GET_T2: Not enough memory.'
                  call xcc_errex
               end if
               call getlst(dcore(ndx_t2(irrep,2)),1,j,2,irrep,listnum)
               freecore = itmp
            end if
         end do
         itmp = freecore - ndx_t2(1,2)
c         dtmp = ddot(itmp,dcore(ndx_t2(1,2)),1,
c     &                    dcore(ndx_t2(1,2)),1)
         dtmp = dnormsqr(itmp,dcore(ndx_t2(1,2)),1)
         t_sqr(2) = t_sqr(2) + dtmp
      end if
#ifdef _DEBUG_XCC
      print *, '<0|adj(T2BB)*T2BB|0> = ', dtmp
#endif

c   o T2(Ij,Ab); Ab, Ij
      listnum = h_T2_IjAb_Ab_Ij
      do irrep = 1, nirrep
         ndx_t2(irrep,3) = freecore
         i = irpdpd(irrep,h_DistType_Ab)
         j = irpdpd(irrep,h_DistType_Ij)
         ioff = i * j
         if (ioff.ne.0) then
            itmp = freecore + ioff
            if (itmp.gt.(1+dcoresize)) then
               print *, '@GET_T2: Not enough memory.'
               call xcc_errex
            end if
            call getlst(dcore(ndx_t2(irrep,3)),1,j,2,irrep,listnum)
            freecore = itmp
         end if
      end do
      itmp = freecore - ndx_t2(1,3)
c      dtmp = ddot(itmp,dcore(ndx_t2(1,3)),1,
c     &                 dcore(ndx_t2(1,3)),1)
      dtmp = dnormsqr(itmp,dcore(ndx_t2(1,3)),1)
      t_sqr(2) = t_sqr(2) + dtmp
#ifdef _DEBUG_XCC
      print *, '<0|adj(T2AB)*T2AB|0> = ', dtmp
#endif /* _DEBUG_XCC */

      overlap = overlap + t_sqr(2)
      print *, '<0|adj(T2)*T2|0> = ', t_sqr(2)

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

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

