
c This routine loads the <AB||IJ> and <Ab|Ij> integrals and contracts them
c with 1/eps(IJ,AB) and 1/eps(Ij,Ab) respectively. These are the first
c approximate T2 amplitudes and are needed for factorized quadruples under CC.

c INPUT
c    double precision dcore(*) : address space for xcc.com

#include "disttype.h"
#include "listnum.h"

c#define _PRINT_AMPLITUDES

      subroutine get_t2p1(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c INTERNAL VARIABLES
      integer listnum, i, j
      integer irrep, ioff, itmp
      double precision dtmp, mbpt2

c PARAMETERS
      double precision zero
      parameter (zero=0.0d0)

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"
#include "method.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      mbpt2 = zero

c   o T2(1)(IJ,AB); A<B, I<J
      listnum = h_HH_ABIJ_AltB_IltJ
      dtmp = zero
      do irrep = 1, nirrep
         ndx_t2p1(irrep,1) = freecore
         i = irpdpd(irrep,h_DistType_AltB)
         j = irpdpd(irrep,h_DistType_IltJ)
         ioff = i * j
         if (ioff.ne.0) then
            itmp = freecore + ioff
            if (itmp.gt.(1+dcoresize)) then
               write(*,*) '@GET_T2P1: Not enough memory.'
               call xcc_errex
            end if
            call getlst(dcore(ndx_t2p1(irrep,1)),1,j,2,irrep,listnum)
            freecore = itmp
            call cont_eps_alta_ilti(
     &              dcore(ndx_t2p1(irrep,1)),
     &              irrep,vrt(1,1),dcore,ndx_eval(1,2,1),i,
     &              irrep,pop(1,1),dcore,ndx_eval(1,1,1),
     &              dtmp
     &           )
         end if
      end do
      mbpt2 = mbpt2 + dtmp
#ifdef _DEBUG_XCC
      write(*,*) 'MBPT(2): E[2]AA = ', dtmp
#endif /* _DEBUG_XCC */

c   o T2(1)(ij,ab); a<b, i<j
      if (iuhf.eq.0) then
         mbpt2 = mbpt2 + dtmp
         do irrep = 1, nirrep
            ndx_t2p1(irrep,2) = ndx_t2p1(irrep,1)
         end do
      else
         listnum = h_HH_abij_altb_iltj
         dtmp = zero
         do irrep = 1, nirrep
            ndx_t2p1(irrep,2) = freecore
            i = irpdpd(irrep,h_DistType_altb)
            j = irpdpd(irrep,h_DistType_iltj)
            ioff = i * j
            if (ioff.ne.0) then
               itmp = freecore + ioff
               if (itmp.gt.(1+dcoresize)) then
                  write(*,*) '@GET_T2P1: Not enough memory.'
                  call xcc_errex
               end if
               call getlst(dcore(ndx_t2p1(irrep,2)),1,j,2,irrep,listnum)
               freecore = itmp
               call cont_eps_alta_ilti(
     &                 dcore(ndx_t2p1(irrep,2)),
     &                 irrep,vrt(1,2),dcore,ndx_eval(1,2,2),i,
     &                 irrep,pop(1,2),dcore,ndx_eval(1,1,2),
     &                 dtmp
     &              )
            end if
         end do
         mbpt2 = mbpt2 + dtmp
      end if
#ifdef _DEBUG_XCC
      write(*,*) 'MBPT(2): E[2]BB = ', dtmp
#endif

c   o T2(Ij,Ab); Ab, Ij
      listnum = h_GS_AbIj_Ab_Ij
      dtmp = zero
      do irrep = 1, nirrep
         ndx_t2p1(irrep,3) = freecore
         i = irpdpd(irrep,h_DistType_Ab)
         j = irpdpd(irrep,h_DistType_Ij)
         ioff = i * j
         if (ioff.ne.0) then
            itmp = freecore + ioff
            if (itmp.gt.(1+dcoresize)) then
               write(*,*) '@GET_T2P1: Not enough memory.'
               call xcc_errex
            end if
            call getlst(dcore(ndx_t2p1(irrep,3)),1,j,2,irrep,listnum)
            freecore = itmp
            call cont_eps_ab_ij(
     &              dcore(ndx_t2p1(irrep,3)),
     &              irrep,vrt(1,1),dcore,ndx_eval(1,2,1),
     &                    vrt(1,2),dcore,ndx_eval(1,2,2),i,
     &              irrep,pop(1,1),dcore,ndx_eval(1,1,1),
     &                    pop(1,2),dcore,ndx_eval(1,1,2),
     &              dtmp
     &           )
         end if
      end do
      mbpt2 = mbpt2 + dtmp
#ifdef _DEBUG_XCC
      write(*,*) 'MBPT(2): E[2]AB = ', dtmp
#endif /* _DEBUG_XCC */

      write(*,*) 'MBPT(2): E[2] = ', mbpt2

#ifdef _PRINT_AMPLITUDES
      do irrep = 1, nirrep
         i = irpdpd(irrep,h_DistType_AltB)
         j = irpdpd(irrep,h_DistType_IltJ)
         write(*,*) 'There are ',i*j,
     &      ' T2(1)AA amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_t2p1(irrep,1)),i,j,i)
      end do
      if (iuhf.eq.1) then
         do irrep = 1, nirrep
            i = irpdpd(irrep,h_DistType_altb)
            j = irpdpd(irrep,h_DistType_iltj)
            write(*,*) 'There are ',i*j,
     &         ' T2(1)BB amplitudes in irrep ',irrep
            call dmat_prt(dcore(ndx_t2p1(irrep,2)),i,j,i)
         end do
      end if
      do irrep = 1, nirrep
         i = irpdpd(irrep,h_DistType_Ab)
         j = irpdpd(irrep,h_DistType_Ij)
         write(*,*) 'There are ',i*j,
     &      ' T2(1)AB amplitudes in irrep ',irrep
         call dmat_prt(dcore(ndx_t2p1(irrep,3)),i,j,i)
      end do
#endif /* _PRINT_AMPLITUDES */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

