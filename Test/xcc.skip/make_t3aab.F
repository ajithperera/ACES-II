
c This routine builds the T3AAB and T3BBA amplitudes. If the reference is
c RHF, then all sorts of ugly things must be done including making the
c T3(I<J<K,A<B<C) amplitudes. Pay attention to the "if (iuhf.eq.0) then"
c lines.

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    integer          ispin1   : the primary spin case to build
c    double precision e4taab   : the E[4]TAAB/BBA correlation energy
c    double precision e5staab  : the E[5]TAAB/BBA correlation energy
c    double precision e4taaa   : the E[4]TAAA/BBB correlation energy
c    double precision e5staaa  : the E[5]TAAA/BBB correlation energy

c OUTPUT
c    double precision t3aab_ovrlp : <0|adj(T3AAB)*T3AAB|0> (not zeroed)
c    double precision t3aaa_ovrlp : <0|adj(T3AAA)*T3AAA|0> (not zeroed)

#include "disttype.h"
#include "listnum.h"
#include "xcc.h"

c#define _DEBUG_MAKE_T3AAB
#define _USE_MPI_BATCH
#define _DO_ABCI
#define _DO_T2IJ /* requires _DO_ABCI */
#define _DO_T2KJ /* requires _DO_ABCI */
#define _DO_T2KI /* requires _DO_ABCI */
#define _DO_T2IK /* requires _DO_ABCI */
#define _DO_T2JK /* requires _DO_ABCI */
#define _DO_IJKA
#define _DO_HHIJ /* requires _DO_IJKA */
#define _DO_GSKJ /* requires _DO_IJKA */
#define _DO_GSKI /* requires _DO_IJKA */
#define _DO_GSJK /* requires _DO_IJKA */
#define _DO_GSIK /* requires _DO_IJKA */
#define _DO_T3AAA_RHF
#define _PROCESS_T3AAA
#define _DO_EPS
#define _PROCESS_T3AAB

      subroutine make_t3aab(dcore,ispin1,
     &                      e4taab,e5staab,t3aab_ovrlp,
     &                      e4taaa,e5staaa,t3aaa_ovrlp)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)
      double precision e4taab, e5staab, t3aab_ovrlp
      double precision e4taaa, e5staaa, t3aaa_ovrlp
      integer ispin1

c INTERNAL VARIABLES
      integer ispin2
      integer start_t3, start_scr, dim_scr
      integer list_i1j1k1a1, dist_i1a1, dist_i1ltj1
      integer list_i1j2k1a2, dist_i1a2, dist_i1j2
      integer list_i2j1k2a1, dist_i2a1
      integer list_a1b1c1i1, dist_a1ltb1
      integer list_a1b2c1i2, dist_a1b2
      integer list_a2b1c2i1
      integer start_t2ij, start_t2jk, start_t2ik
      integer add_t2ij,   add_t2jk,   add_t2ik
      integer start_wik, start_wjk, start_wij, start_wkj, start_wki
      integer add_wik,   add_wjk,   add_wij,   add_wkj,   add_wki
      integer i, irp_i, max_i,         to_i
      integer j, irp_j, max_j, from_j, to_j
      integer k, irp_k, max_k
      integer irp_ij, irp_jk, irp_ik, irp_ijk
      logical do_irp_j
      integer irrep, itmp, ioff
      integer ndx_gs21(8), ndx_gs12(8), ndx_hh11(8)
      integer off_i1ltj1, off_a1ltb1
      integer off_i2j1, off_i1j2, off_a1b2, off_a1i1
      integer off_altbltc(37)
      integer off_a1ltb1c2(65), ndx_t3vrt(64), max_a1b1c2, max_a1b1c1
      double precision eps_k
      double precision eps_jk
      double precision eps_ijk
      integer mpi_count, mpi_batch

c PARAMETERS
      double precision one, minus
      parameter (one=1.0d0, minus=-1.0d0)

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"
#include "method.com"
#include "machsp.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
c   o assert ispin1 is 1 or 2
      if ((ispin1.lt.1).or.(2.lt.ispin1)) then
         print *, '@MAKE_T3AAB: Assertion failed.'
         print *, '   ispin1 = ',ispin1
         call xcc_errex
      end if
#endif /* _ASSERT */

c The plan is:
c   o loop over I<J,k
c      - create all W(+1)*T2 diagrams
c   o divide each term by epsilon(IJk,ABc)

c   o set handles for these dependents
      if (ispin1.eq.1) then
         ispin2        = 2
         off_i1ltj1    = h_Offsets_IltJ
         off_a1ltb1    = h_Offsets_AltB
         off_i1j2      = h_Offsets_Ij
         off_i2j1      = h_Offsets_iJ
         off_a1i1      = h_Offsets_AI
         off_a1b2      = h_Offsets_Ab
         list_i1j1k1a1 = h_HH_IJKA_IltJ_KA
         list_i1j2k1a2 = h_GS_IjKa_Ij_Ka
         list_i2j1k2a1 = h_GS_IjAk_Ij_Ak
         list_a1b1c1i1 = h_HH_ABCI_AltB_CI
         list_a1b2c1i2 = h_GS_AbCi_Ab_Ci
         list_a2b1c2i1 = h_GS_AbIc_Ab_Ic
         dist_i1ltj1   = h_DistType_IltJ
         dist_i1a1     = h_DistType_IA
         dist_i1a2     = h_DistType_Ia
         dist_i2a1     = h_DistType_Ai
         dist_a1ltb1   = h_DistType_AltB
      else
         ispin2        = 1
         off_i1ltj1    = h_Offsets_iltj
         off_a1ltb1    = h_Offsets_altb
         off_i1j2      = h_Offsets_iJ
         off_i2j1      = h_Offsets_Ij
         off_a1i1      = h_Offsets_ai
         off_a1b2      = h_Offsets_aB
         list_i1j1k1a1 = h_HH_ijka_iltj_ka
         list_i1j2k1a2 = h_GS_IjAk_Ij_Ak
         list_i2j1k2a1 = h_GS_IjKa_Ij_Ka
         list_a1b1c1i1 = h_HH_abci_altb_ci
         list_a1b2c1i2 = h_GS_AbIc_Ab_Ic
         list_a2b1c2i1 = h_GS_AbCi_Ab_Ci
         dist_i1ltj1   = h_DistType_iltj
         dist_i1a1     = h_DistType_ia
         dist_i1a2     = h_DistType_Ai
         dist_i2a1     = h_DistType_Ia
         dist_a1ltb1   = h_DistType_altb
      end if
      dist_i1j2 = h_DistType_Ij
      dist_a1b2 = h_DistType_Ab

c   o initialize some constants and variables
      start_t3  = ndx_t3(1,2+ispin1)
      start_scr = freecore
      dim_scr   = 1 + dcoresize - freecore

c   o get all the <Ij|Ka> integrals and transpose them to <Ka|Ij>
c     NOTE: If (ispin1.eq.2), the list is actually <Ij|Ak> so each
c           half of the integral must be individually transposed.
c     NOTE: If this is an RHF calculation, then these integrals will
c           be used to make <IJ||KA> as soon as they are read in.
c     start_scr will eventually point to the address after the integrals.
      do irrep = 1, nirrep
         ndx_gs12(irrep) = start_scr
         i = irpdpd(irrep,dist_i1j2)
         j = irpdpd(irrep,dist_i1a2)
         ioff = i * j
         if (ioff.ne.0) then
            start_scr = start_scr + ioff
            dim_scr   = 1 + dcoresize - start_scr
            if (dim_scr.lt.1) then
               print *, '@MAKE_T3AAB: Not enough memory.'
               call xcc_errex
            end if
            call getlst(dcore(ndx_gs12(irrep)),1,j,2,irrep,
     &                  list_i1j2k1a2)
            if (iuhf.eq.0) then
c            o ioff is now used to measure the number of <IJ||KA> integrals
               ndx_hh11(irrep) = start_scr
               k = irpdpd(irrep,dist_i1ltj1)
               ioff = k * j
               if (ioff.ne.0) then
                  start_scr = start_scr + ioff
                  dim_scr   = 1 + dcoresize - start_scr
                  if (dim_scr.lt.1) then
                     print *, '@MAKE_T3AAB: Not enough memory.'
                     call xcc_errex
                  end if
                  call antisym_pp_pltp(dcore(ndx_gs12(irrep)),i,
     &                                 dcore(ndx_hh11(irrep)),
     &                                 pop(1,1),j,
     &                                 irrep,off_pq(1,off_i1j2))
                  if ((k.gt.1).and.(j.gt.1)) then
                     call dmat_xpose(dcore(ndx_hh11(irrep)),k,j)
                  end if
               end if
            end if
            if (ispin1.eq.2) then
               call transp_dist(dcore(ndx_gs12(irrep)),
     &                          pop(1,1),pop(1,2),
     &                          irrep,j,.false.)
            end if
            if ((i.gt.1).and.(j.gt.1)) then
               call dmat_xpose(dcore(ndx_gs12(irrep)),i,j)
            end if
            if (ispin1.eq.2) then
               call transp_dist(dcore(ndx_gs12(irrep)),
     &                          vrt(1,1),pop(1,2),
     &                          irrep,i,.false.)
            end if
c        else if (ioff.eq.0) then
         else
            if (iuhf.eq.0) then
c            o make sure ndx_hh11(irrep) is non-zero
               ndx_hh11(irrep) = start_scr
            end if
         end if
      end do

c   o get all the extra UHF integrals
      if (iuhf.ne.0) then

c      o get all the <iJ|kA> integrals and transpose them to <kA|iJ>
c        NOTE: If (ispin1.eq.1), the list is actually <Ij|Ak> so each
c              half of the integral must be individually transposed.
c              If, however, this is an RHF calculation, then the <Ka|Ij>
c              list in memory will be used.
c        start_scr will eventually point to the address after the integrals.
         do irrep = 1, nirrep
            ndx_gs21(irrep) = start_scr
            i = irpdpd(irrep,dist_i1j2)
            j = irpdpd(irrep,dist_i2a1)
            ioff = i * j
            if (ioff.ne.0) then
               start_scr = start_scr + ioff
               dim_scr   = 1 + dcoresize - start_scr
               if (dim_scr.lt.1) then
                  print *, '@MAKE_T3AAB: Not enough memory.'
                  call xcc_errex
               end if
               call getlst(dcore(ndx_gs21(irrep)),1,j,2,irrep,
     &                     list_i2j1k2a1)
               if (ispin1.eq.1) then
                  call transp_dist(dcore(ndx_gs21(irrep)),
     &                             pop(1,1),pop(1,2),
     &                             irrep,j,.false.)
               end if
               if ((i.gt.1).and.(j.gt.1)) then
                  call dmat_xpose(dcore(ndx_gs21(irrep)),i,j)
               end if
               if (ispin1.eq.1) then
                  call transp_dist(dcore(ndx_gs21(irrep)),
     &                             vrt(1,1),pop(1,2),
     &                             irrep,i,.false.)
               end if
            end if
         end do

c      o get all the <IJ||KA> integrals and transpose them to <KA||IJ>
c        start_scr will eventually point to the address after the integrals.
         do irrep = 1, nirrep
            ndx_hh11(irrep) = start_scr
            i = irpdpd(irrep,dist_i1ltj1)
            j = irpdpd(irrep,dist_i1a1)
            ioff = i * j
            if (ioff.ne.0) then
               start_scr = start_scr + ioff
               dim_scr   = 1 + dcoresize - start_scr
               if (dim_scr.lt.1) then
                  print *, '@MAKE_T3AAB: Not enough memory.'
                  call xcc_errex
               end if
               call getlst(dcore(ndx_hh11(irrep)),1,j,2,irrep,
     &                     list_i1j1k1a1)
               if ((i.gt.1).and.(j.gt.1)) then
                  call dmat_xpose(dcore(ndx_hh11(irrep)),i,j)
               end if
            end if
         end do

c     else if (iuhf.eq.0) then
      else
         do irrep = 1, nirrep
            ndx_gs21(irrep) = ndx_gs12(irrep)
         end do
c      o The <IJ||KA> integrals were already created and transposed.
c     end if (iuhf.ne.0)
      end if

c ----------------------------------------------------------------------

      mpi_count = 0
      mpi_batch = (224+CNTOT-1)/CNTOT

      do irp_ijk = 1, nirrep

c      o off_altbltc and off_a1ltb1c2 contain the offsets of virtuals
c        needed to address irp_a, irp_b, and irp_c whose direct product
c        is irp_ijk.
c        NOTE: THESE ELEMENTS DO NOT BEHAVE LIKE THE ONES IN off_pltp.
         if (iuhf.eq.0) then
            call offset_cum_pltpltp(irp_ijk,vrt(1,1),off_altbltc)
            max_a1b1c1 = off_altbltc(1+nirrep*(nirrep+1)/2)
         end if
         call offset_cum_pltpq(irp_ijk,vrt(1,ispin1),vrt(1,ispin2),
     &                         off_a1ltb1c2)
         max_a1b1c2 = off_a1ltb1c2(nirrep*nirrep+1)
      if (max_a1b1c2.ne.0) then

c   o loop over irreps of occupieds in T3(IJk,ABc)
      do irp_k = 1, nirrep
         max_k = pop(irp_k,ispin2)
      if (max_k.ne.0) then
      do irp_j = 1, nirrep
         irp_jk = dirprd(irp_j,irp_k)
         irp_i  = dirprd(irp_jk,irp_ijk)
         if (irp_i.lt.irp_j) then
            max_i = pop(irp_i,ispin1)
            max_j = pop(irp_j,ispin1)
            do_irp_j = ((max_j.ne.0).and.(max_i.ne.0))
         else
            if (irp_i.gt.irp_j) then
               do_irp_j = .false.
            else
               max_i = pop(irp_i,ispin1)
               max_j = max_i
               do_irp_j = (max_j.gt.1)
            end if
         end if
      if (do_irp_j) then

c      o evaluate symmetries once
         irp_ij = dirprd(irp_j,irp_i)
         irp_ik = dirprd(irp_k,irp_i)

c      o set orbital boundary conditions (as offsets)
         to_i = max_i-1
         to_j = max_j-1
         if (irp_i.eq.irp_j) then
            from_j = 1
         else
            from_j = 0
         end if

c      o find the T2(IJ,*) amplitudes and <*||IJ> integrals
         itmp = rshift(irp_j*(irp_j-1),1) + irp_i
         ioff = off_pltp(itmp,off_i1ltj1)
         start_t2ij =   ndx_t2(irp_ij,ispin1)
     &                + ( ioff * irpdpd(irp_ij,dist_a1ltb1) )
         start_wij  =   ndx_hh11(irp_ij)
     &                + ( ioff * irpdpd(irp_ij,dist_i1a1) )

c      o find the <*|Ik> and <*|Jk> integrals
         itmp = nirrep*(irp_k-1)
         i = itmp + irp_i
         j = itmp + irp_j
         i = off_pq(i,off_i1j2)
         j = off_pq(j,off_i1j2) + from_j
         start_wik =   ndx_gs12(irp_ik)
     &               + ( i * irpdpd(irp_ik,dist_i1a2) )
         start_wjk =   ndx_gs12(irp_jk)
     &               + ( j * irpdpd(irp_jk,dist_i1a2) )

c      o find the <*|kI> and <*|kJ> integrals
         i = nirrep*(irp_i-1) + irp_k
         j = nirrep*(irp_j-1) + irp_k
         i = off_pq(i,off_i2j1)
         j = off_pq(j,off_i2j1) + ( from_j * max_k )
         start_wki =   ndx_gs21(irp_ik)
     &               + ( i * irpdpd(irp_ik,dist_i2a1) )
         start_wkj =   ndx_gs21(irp_jk)
     &               + ( j * irpdpd(irp_jk,dist_i2a1) )

c      o find the T2(Ik,*) and T2(Jk,*) amplitudes
c        All these have to do are point to the correct Ik/Jk lists.
c        The subroutine(s) will figure out how the amplitudes are
c        actually stored and resort/scale them accordingly.
         if (ispin1.eq.1) then
            itmp = nirrep*(irp_k-1)
            i = itmp + irp_i
            j = itmp + irp_j
            i = off_pq(i,h_Offsets_Ij)
            j = off_pq(j,h_Offsets_Ij) + from_j
         else
            i = nirrep*(irp_i-1) + irp_k
            j = nirrep*(irp_j-1) + irp_k
            i = off_pq(i,h_Offsets_Ij)
            j = off_pq(j,h_Offsets_Ij) + ( from_j * max_k )
         end if
         start_t2ik =   ndx_t2(irp_ik,3)
     &                + ( i * irpdpd(irp_ik,dist_a1b2) )
         start_t2jk =   ndx_t2(irp_jk,3)
     &                + ( j * irpdpd(irp_jk,dist_a1b2) )

c      o loop over orbitals I<J,k
            add_wjk  = start_wjk
         do k = 0, max_k-1
            eps_k = dcore(ndx_eval(irp_k,1,ispin2) + k)
            add_t2ij = start_t2ij
            add_t2jk = start_t2jk
            add_wij  = start_wij
            add_wkj  = start_wkj
         do j = from_j, to_j
            eps_jk = eps_k + dcore(ndx_eval(irp_j,1,ispin1) + j)
            add_t2ik = start_t2ik
            add_wik  = start_wik
            add_wki  = start_wki
            if (irp_i.eq.irp_j) to_i = j-1
         do i = 0, to_i

      mpi_count = mpi_count + 1
#ifdef _USE_MPI_BATCH
      if (((CNID*mpi_batch).lt.mpi_count).and.
     &                        (mpi_count.le.((CNID+1)*mpi_batch))) then
#else
      if (mod(mpi_count,CNTOT).eq.CNID) then
#endif

c ----------------------------------------------------------------------

      eps_ijk = eps_jk + dcore(ndx_eval(irp_i,1,ispin1) + i)

c   o build an address list for dependents
      do itmp = 1, nirrep*nirrep
         ndx_t3vrt(itmp) = start_t3 + off_a1ltb1c2(itmp)
      end do

c   o initialize the T3 amplitudes
      do itmp = 0, max_a1b1c2-1
         dcore(start_t3+itmp) = 0.0d0
      end do

#ifdef _DO_ABCI
#ifdef _DO_T2IJ
c   o   <Pq|Ek> * T2(IJ,ER)
      call cont_gs_a1b2c1i2_t2aa(
     &        dcore,  dcore(start_scr), dim_scr,
     &        ispin1, ispin2,
     &        list_a1b2c1i2, dist_a1b2,
     &        add_t2ij, off_a1ltb1,
     &        i, irp_i, max_i, irp_ij,
     &        j, irp_j, max_j, irp_jk,
     &        k, irp_k, max_k, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_T2KJ
c   o   <pQ|eI> * T2(kJ,eR)
      call cont_gs_a2b1c2i1_t2ab(
     &        dcore,  dcore(start_scr), dim_scr,
     &        ispin1, ispin2,
     &        list_a2b1c2i1, dist_a1b2,
     &        i, irp_i, max_i,
     &        add_t2jk, off_a1b2,
     &        j, irp_j, max_j, irp_jk,
     &        k, irp_k, max_k, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_T2KI
c   o - <pQ|eJ> * T2(kI,eR)
      call cont_gs_a2b1c2i1_t2ab(
     &        dcore,  dcore(start_scr), dim_scr,
     &        ispin1, ispin2,
     &        list_a2b1c2i1, dist_a1b2,
     &        j, irp_j, max_j,
     &        add_t2ik, off_a1b2,
     &        i, irp_i, max_i, irp_ik,
     &        k, irp_k, max_k, irp_ijk,
     &        minus, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_T2JK
c   o   <PQ||EI> * T2(Jk,Er)
      call cont_hh_abci_t2ab(
     &        dcore,  dcore(start_scr), dim_scr,
     &        ispin1, ispin2,
     &        list_a1b1c1i1, dist_a1ltb1, off_a1i1,
     &        i, irp_i, max_i,
     &        add_t2jk, off_a1b2,
     &        j, irp_j, max_j, irp_jk,
     &        k, irp_k, max_k, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_T2IK
c   o - <PQ||EJ> * T2(Ik,Er)
      call cont_hh_abci_t2ab(
     &        dcore,  dcore(start_scr), dim_scr,
     &        ispin1, ispin2,
     &        list_a1b1c1i1, dist_a1ltb1, off_a1i1,
     &        j, irp_j, max_j,
     &        add_t2ik, off_a1b2,
     &        i, irp_i, max_i, irp_ik,
     &        k, irp_k, max_k, irp_ijk,
     &        minus, start_t3, ndx_t3vrt
     &     )
#endif
#endif /* _DO_ABCI */

#ifdef _DO_IJKA
#ifdef _DO_HHIJ
c   o   T2(Mk,Pq) * <MR||IJ>
      call cont_hh_iajk_t2ab(
     &        dcore,    dcore(start_scr), dim_scr,
     &        ispin1,   ispin2,
     &        off_i1j2, dist_a1b2,
     &        k, irp_k, max_k,
     &        add_wij,
     &        i, irp_i, max_i, irp_ij,
     &        j, irp_j, max_j, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_GSKJ
c   o   T2(mI,pQ) * <mR|kJ>
      call cont_gs_i2a1j2k1_t2ab(
     &        dcore,  dcore(start_scr), dim_scr,
     &        ispin1, ispin2,
     &        off_i1j2, dist_a1b2,
     &        i, irp_i, max_i,
     &        add_wkj,
     &        k, irp_k, max_k, irp_jk,
     &        j, irp_j, max_j, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_GSKI
c   o - T2(mJ,pQ) * <mR|kI>
      call cont_gs_i2a1j2k1_t2ab(
     &        dcore,  dcore(start_scr), dim_scr,
     &        ispin1, ispin2,
     &        off_i1j2, dist_a1b2,
     &        j, irp_j, max_j,
     &        add_wki,
     &        k, irp_k, max_k, irp_ik,
     &        i, irp_i, max_i, irp_ijk,
     &        minus, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_GSJK
c   o   T2(MI,PQ) * <Mr|Jk>
      call cont_gs_i1a2j1k2_t2aa(
     &        dcore,  dcore(start_scr), dim_scr,
     &        ispin1, ispin2,
     &        off_i1ltj1, dist_a1ltb1,
     &        add_wjk,
     &        i, irp_i, max_i, irp_ij,
     &        j, irp_j, max_j, irp_jk,
     &        k, irp_k, max_k, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_GSIK
c   o - T2(MJ,PQ) * <Mr|Ik>
      call cont_gs_i1a2j1k2_t2aa(
     &        dcore,  dcore(start_scr), dim_scr,
     &        ispin1, ispin2,
     &        off_i1ltj1, dist_a1ltb1,
     &        add_wik,
     &        j, irp_j, max_j, irp_ij,
     &        i, irp_i, max_i, irp_ik,
     &        k, irp_k, max_k, irp_ijk,
     &        minus, start_t3, ndx_t3vrt
     &     )
#endif
#endif /* _DO_IJKA */

#ifdef _DO_T3AAA_RHF
      if ((iuhf.eq.0).and.
     &    ((irp_j.lt.irp_k).or.
     &     ((irp_j.eq.irp_k).and.(j.lt.k))
     &    )
     &   ) then
c      o make the T3AAA amplitudes
         if (max_a1b1c1.ne.0) then
            call make_t3aaa_rhf(
     &              dcore,ndx_t3vrt,
     &              dcore(start_scr),off_altbltc,
     &              dcore(start_scr+max_a1b1c1),dim_scr-max_a1b1c1,
     &              i,irp_i,max_i,irp_ij,
     &              j,irp_j,max_j,irp_jk,
     &              k,irp_k,max_k,irp_ik,irp_ijk,
     &              eps_ijk,e4taaa,t3aaa_ovrlp
     &           )
#ifdef _PROCESS_T3AAA
            call process_t3aaa(
     &              dcore,
     &              dcore(start_scr+max_a1b1c1),dim_scr-max_a1b1c1,
     &              ispin1,
     &              dcore(start_scr),max_a1b1c1,off_altbltc,
     &              i,irp_i,max_i,irp_ij,add_t2ij,
     &              j,irp_j,max_j,irp_jk,add_t2jk,
     &              k,irp_k,max_k,irp_ik,add_t2ik,irp_ijk,
     &              e5staaa
     &           )
#endif /* _PROCESS_T3AAA */
         end if
      end if
#endif /* _DO_T3AAA_RHF */

#ifdef _DO_EPS
c   o divide the pre-amplitudes with epsilon(IJk,ABc) and return
c     E[4]T and the T3 overlap contribution
      call cont_eps_t3aab(dcore,dcore(start_t3),
     &                    eps_ijk,irp_ijk,
     &                    ispin1,ispin2,
     &                    e4taab,t3aab_ovrlp)
#endif /* _DO_EPS */

#ifdef _PROCESS_T3AAB
      call process_t3aab(
     &        dcore,dcore(start_scr+max_a1b1c2),dim_scr-max_a1b1c2,
     &        ispin1,
     &        dcore(start_t3),max_a1b1c2,off_a1ltb1c2,
     &        i,irp_i,max_i,irp_ij,
     &        j,irp_j,max_j,irp_jk,
     &        k,irp_k,max_k,irp_ik,irp_ijk,
     &        e5staab
     &     )
#endif /* _PROCESS_T3AAB */

#ifndef ONLY_ONE_T3_DIST
c   o point to the next T3(IJk,*) distribution
      start_t3 = start_t3 + max_a1b1c2
#endif

c ----------------------------------------------------------------------

c     end if (mod(mpi_count,20).eq.CNID)
      end if

c      o increment and close the loop over orbitals I<J,k
            add_t2ij = add_t2ij + irpdpd(irp_ij,dist_a1ltb1)
            add_wki  = add_wki  + ( irpdpd(irp_ik,dist_i2a1) * max_k )
            if (ispin1.eq.1) then
               add_t2ik = add_t2ik + irpdpd(irp_ik,dist_a1b2)
            else
               add_t2ik =   add_t2ik
     &                    + ( irpdpd(irp_ik,dist_a1b2) * max_k )
            end if
            add_wij = add_wij + irpdpd(irp_ij,dist_i1a1)
            add_wik = add_wik + irpdpd(irp_ik,dist_i1a2)
c        end do i = 0, to_i
         end do
            add_wkj = add_wkj + ( irpdpd(irp_jk,dist_i2a1) * max_k )
            if (ispin1.eq.1) then
               add_t2jk = add_t2jk + irpdpd(irp_jk,dist_a1b2)
            else
               add_t2jk =   add_t2jk
     &                    + ( irpdpd(irp_jk,dist_a1b2) * max_k )
            end if
            add_wjk = add_wjk + irpdpd(irp_jk,dist_i1a2)
c        end do j = from_j, to_j
         end do
            start_wki = start_wki + irpdpd(irp_ik,dist_i2a1)
            start_wkj = start_wkj + irpdpd(irp_jk,dist_i2a1)
            if (ispin1.eq.1) then
               start_t2ik =   start_t2ik
     &                      + ( irpdpd(irp_ik,dist_a1b2) * max_i )
               start_t2jk =   start_t2jk
     &                      + ( irpdpd(irp_jk,dist_a1b2) * max_j  )
            else
               start_t2ik = start_t2ik + irpdpd(irp_ik,dist_a1b2)
               start_t2jk = start_t2jk + irpdpd(irp_jk,dist_a1b2)
            end if
            add_wjk   = add_wjk + ( irpdpd(irp_jk,dist_i1a2) * from_j )
            start_wik =   start_wik
     &                  + ( irpdpd(irp_ik,dist_i1a2) * max_i )
c        end do k = 0, max_k-1
         end do

c   o close the loop over irreps of occupieds in T3(IJk,ABc)
c     end if (do_irp_j)
      end if
c     end do irp_j = 1, nirrep
      end do
c     end if (max_k.ne.0)
      end if
c     end do irp_k = 1, nirrep
      end do

c     end if (max_a1b1c2.ne.0)
      end if
c     end do irp_ijk = 1, nirrep
      end do

      if (ispin1.eq.1) then
         print *, 'There are ',mpi_count,' AAB triples.'
      else
         print *, 'There are ',mpi_count,' BBA triples.'
      end if

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

