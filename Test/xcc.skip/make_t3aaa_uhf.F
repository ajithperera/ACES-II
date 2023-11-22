
c This routine builds the T3AAA and T3BBB amplitudes for UHF calculations.

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    integer          ispin    : the spin case to build
c    double precision e4t      : the E[4]TAAA/BBB correlation energy
c    double precision e5st     : the E[5]TAAA/BBB correlation energy

c OUTPUT
c    double precision t3_ovrlp : <0|adj(T3AAA)*T3AAA|0> (not zeroed)

#include "disttype.h"
#include "listnum.h"
#include "xcc.h"

c#define _DEBUG_MAKE_T3AAA_UHF
#define _USE_MPI_BATCH
#define _DO_ABCI
#define _DO_T2_JK /* requires _DO_ABCI */
#define _DO_T2_IK /* requires _DO_ABCI */
#define _DO_T2_IJ /* requires _DO_ABCI */
#define _DO_IJKA
#define _DO_HH_JK /* requires _DO_IJKA */
#define _DO_HH_IK /* requires _DO_IJKA */
#define _DO_HH_IJ /* requires _DO_IJKA */
#define _DO_EPS
#define _PROCESS_T3AAA

      subroutine make_t3aaa_uhf(dcore,ispin,e4t,e5st,t3_ovrlp)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), e4t, e5st, t3_ovrlp
      integer ispin

c INTERNAL VARIABLES
      integer list_abci, list_ijka
      integer dist_altb, dist_iltj, dist_ia
      integer off_iltj, off_altb, off_ai
      integer i, irp_i, max_i,         to_i
      integer j, irp_j, max_j, from_j, to_j
      integer k, irp_k, max_k, from_k
      integer irp_ij, irp_ik, irp_jk, irp_ijk
      logical do_irp_j
      integer start_scr,  dim_scr
      integer start_wij,  add_wij
      integer start_wik,  add_wik
      integer start_wjk,  add_wjk, ndx_w(8)
      integer start_t2ij, add_t2ij
      integer start_t2ik, add_t2ik
      integer start_t2jk, add_t2jk
      integer start_t3,   off_altbltc(37), ndx_t3vrt(36), dim_upper_tri
      integer max_abc
      integer irrep, itmp, ioff
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
#include "machsp.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
c   o assert ispin is 1 or 2
      if ((ispin.lt.1).or.(2.lt.ispin)) then
         print *, '@MAKE_T3AAA_UHF: Assertion failed.'
         print *, '   ispin = ',ispin
         call xcc_errex
      end if
#endif /* _ASSERT */

c The plan is:
c   o loop over I<J<K
c      - create all W(+1)*T2 diagrams
c   o divide each term by epsilon(IJK,ABC)

c   o set handles for these dependents
      if (ispin.eq.1) then
         list_ijka = h_HH_IJKA_IltJ_KA
         list_abci = h_HH_ABCI_AltB_CI
         dist_altb = h_DistType_AltB
         dist_iltj = h_DistType_IltJ
         dist_ia   = h_DistType_IA
         off_iltj  = h_Offsets_IltJ
         off_altb  = h_Offsets_AltB
         off_ai    = h_Offsets_AI
      else
         list_ijka = h_HH_ijka_iltj_ka
         list_abci = h_HH_abci_altb_ci
         dist_altb = h_DistType_altb
         dist_iltj = h_DistType_iltj
         dist_ia   = h_DistType_ia
         off_iltj  = h_Offsets_iltj
         off_altb  = h_Offsets_altb
         off_ai    = h_Offsets_ai
      end if

c   o initialize some constants and variables
      start_t3  = ndx_t3(1,ispin)
      start_scr = freecore
      dim_scr   = 1 + dcoresize - freecore
      dim_upper_tri = rshift(nirrep*(nirrep+1),1)

c   o get all the <IJ||KA> integrals and transpose them to <KA||IJ>
c     start_scr will eventually point to the address after the OOOV
c     integrals.
      do irrep = 1, nirrep
         ndx_w(irrep) = start_scr
         i = irpdpd(irrep,dist_iltj)
         j = irpdpd(irrep,dist_ia)
         ioff = i * j
         if (ioff.ne.0) then
            start_scr = start_scr + ioff
            dim_scr   = 1 + dcoresize - start_scr
            if (dim_scr.lt.1) then
               print *, '@MAKE_T3AAA_UHF: Not enough memory.'
               call xcc_errex
            end if
            call getlst(dcore(ndx_w(irrep)),1,j,2,irrep,list_ijka)
            if ((i.gt.1).and.(j.gt.1)) then
               call dmat_xpose(dcore(ndx_w(irrep)),i,j)
            end if
         end if
      end do

c ----------------------------------------------------------------------

      mpi_count = 0
      mpi_batch = (56+CNTOT-1)/CNTOT

      do irp_ijk = 1, nirrep
c      o off_altbltc contains the offsets of virtuals needed to address
c        irp_a, irp_b, and irp_c whose direct product is irp_ijk.
c        NOTE: THESE ELEMENTS DO NOT BEHAVE LIKE THE ONES IN off_pltp.
         call offset_cum_pltpltp(irp_ijk,vrt(1,ispin),off_altbltc)
         max_abc = off_altbltc(dim_upper_tri+1)
      if (max_abc.ne.0) then

c   o loop over irreps of occupieds in T3(IJK,ABC)
      do irp_k = 1, nirrep
         max_k = pop(irp_k,ispin)
      if (max_k.ne.0) then
      do irp_j = 1, irp_k
         irp_jk = dirprd(irp_j,irp_k)
         irp_i  = dirprd(irp_jk,irp_ijk)
         if (irp_i.lt.irp_j) then
            max_i = pop(irp_i,ispin)
            max_j = pop(irp_j,ispin)
            if (irp_j.eq.irp_k) then
               do_irp_j = ((max_j.ne.1).and.(max_i.ne.0))
            else
               do_irp_j = ((max_j.ne.0).and.(max_i.ne.0))
            end if
         else
            if (irp_i.gt.irp_j) then
               do_irp_j = .false.
            else
               max_i = pop(irp_i,ispin)
               max_j = max_i
               if (irp_j.eq.irp_k) then
                  do_irp_j = (max_j.gt.2)
               else
                  do_irp_j = (max_j.gt.1)
               end if
            end if
         end if
      if (do_irp_j) then

c      o evaluate symmetries once
         irp_ij = dirprd(irp_j,irp_i)
         irp_ik = dirprd(irp_k,irp_i)

c      o set orbital boundary conditions (as offsets)
         to_i = max_i-1
         to_j = max_j-1
         if (irp_i.ne.irp_j) then
            from_j = 0
         else
            from_j = 1
         end if
         if (irp_j.ne.irp_k) then
               from_k = 0
         else
            if (irp_i.ne.irp_j) then
               from_k = 1
            else
               from_k = 2
            end if
         end if

c      o find the T2(IJ,*) amplitudes and <*||IJ> integrals
         itmp = rshift(irp_j*(irp_j-1),1) + irp_i
         ioff = off_pltp(itmp,off_iltj)
         start_t2ij =   ndx_t2(irp_ij,ispin)
     &                + ( ioff * irpdpd(irp_ij,dist_altb) )
         start_wij  =   ndx_w(irp_ij)
     &                + ( ioff * irpdpd(irp_ij,dist_ia) )

c      o find the T2(IK,*) amplitudes and <*||IK> integrals
c        (( these addresses will shift after each k iteration ))
         itmp = rshift(irp_k*(irp_k-1),1) + irp_i
         ioff = off_pltp(itmp,off_iltj)
         if (irp_i.eq.irp_k) then
            ioff = ioff + rshift(from_k*(from_k-1),1)
         else
            ioff = ioff + ( from_k * max_i )
         end if
         start_t2ik =   ndx_t2(irp_ik,ispin)
     &                + ( ioff * irpdpd(irp_ik,dist_altb) )
         start_wik  =   ndx_w(irp_ik)
     &                + ( ioff * irpdpd(irp_ik,dist_ia) )

c      o find the T2(JK,*) amplitudes and <*||JK> integrals
         itmp = rshift(irp_k*(irp_k-1),1) + irp_j
         ioff = off_pltp(itmp,off_iltj)
         if (irp_j.eq.irp_k) then
            ioff = ioff + rshift(from_k*(from_k-1),1) + from_j
         else
            ioff = ioff + ( from_k * max_j ) + from_j
         end if
         start_t2jk =   ndx_t2(irp_jk,ispin)
     &                + ( ioff * irpdpd(irp_jk,dist_altb) )
         start_wjk  =   ndx_w(irp_jk)
     &                + ( ioff * irpdpd(irp_jk,dist_ia) )

c      o loop over orbitals I<J<K
            add_t2jk = start_t2jk
            add_wjk  = start_wjk
         do k = from_k, max_k-1
            eps_k = dcore(ndx_eval(irp_k,1,ispin) + k)
            add_t2ij = start_t2ij
            add_wij  = start_wij
            if (irp_j.eq.irp_k) to_j = k-1
         do j = from_j, to_j
            eps_jk = eps_k + dcore(ndx_eval(irp_j,1,ispin) + j)
            add_t2ik = start_t2ik
            add_wik  = start_wik
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

#ifdef _DEBUG_MAKE_T3AAA_UHF
      print *, '@',__FILE__,': T3( ',
     &         1+i,'[',irp_i,'] ',
     &         1+j,'[',irp_j,'] ',
     &         1+k,'[',irp_k,'] , * )'
#endif

c   o build an address list for dependents
      do itmp = 1, dim_upper_tri
         ndx_t3vrt(itmp) = start_t3 + off_altbltc(itmp)
      end do

c   o initialize the T3 amplitudes
      do itmp = 0, max_abc-1
         dcore(start_t3+itmp) = 0.0d0
      end do

#ifdef _DO_ABCI
#ifdef _DO_T2_JK
c   o   <PQ||EI> * T2(JK,ER)
      call cont_hh_abci_t2aa(
     &        dcore, dcore(start_scr), dim_scr,
     &        ispin,
     &        list_abci, dist_altb, off_ai,
     &        add_t2jk,  off_altb,
     &        i, irp_i,  max_i, irp_ij,
     &        j, irp_j,  max_j, irp_jk,
     &        k, irp_k,  max_k, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_T2_IK
c   o - <PQ||EJ> * T2(IK,ER)
      call cont_hh_abci_t2aa(
     &        dcore, dcore(start_scr), dim_scr,
     &        ispin,
     &        list_abci, dist_altb, off_ai,
     &        add_t2ik,  off_altb,
     &        j, irp_j,  max_j, irp_ij,
     &        i, irp_i,  max_i, irp_ik,
     &        k, irp_k,  max_k, irp_ijk,
     &        minus, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_T2_IJ
c   o   <PQ||EK> * T2(IJ,ER)
      call cont_hh_abci_t2aa(
     &        dcore, dcore(start_scr), dim_scr,
     &        ispin,
     &        list_abci, dist_altb, off_ai,
     &        add_t2ij,  off_altb,
     &        k, irp_k,  max_k, irp_ik,
     &        i, irp_i,  max_i, irp_ij,
     &        j, irp_j,  max_j, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif
#endif /* _DO_ABCI */

#ifdef _DO_IJKA
#ifdef _DO_HH_JK
c   o   T2(MI,PQ) * <MR||JK>
      call cont_hh_iajk_t2aa(
     &        dcore, dcore(start_scr), dim_scr,
     &        ispin,
     &        off_iltj, dist_altb,
     &        i, irp_i, max_i,
     &        add_wjk,
     &        j, irp_j, max_j, irp_jk,
     &        k, irp_k, max_k, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_HH_IK
c   o - T2(MJ,PQ) * <MR||IK>
      call cont_hh_iajk_t2aa(
     &        dcore, dcore(start_scr), dim_scr,
     &        ispin,
     &        off_iltj, dist_altb,
     &        j, irp_j, max_j,
     &        add_wik,
     &        i, irp_i, max_i, irp_ik,
     &        k, irp_k, max_k, irp_ijk,
     &        minus, start_t3, ndx_t3vrt
     &     )
#endif

#ifdef _DO_HH_IJ
c   o   T2(MK,PQ) * <MR||IJ>
      call cont_hh_iajk_t2aa(
     &        dcore, dcore(start_scr), dim_scr,
     &        ispin,
     &        off_iltj, dist_altb,
     &        k, irp_k, max_k,
     &        add_wij,
     &        i, irp_i, max_i, irp_ij,
     &        j, irp_j, max_j, irp_ijk,
     &        one, start_t3, ndx_t3vrt
     &     )
#endif
#endif /* _DO_IJKA */

#ifdef _DO_EPS
c   o divide the pre-amplitudes with epsilon(IJK,ABC) and return
c     E[4]T and the T3 overlap contribution
      eps_ijk = eps_jk + dcore(ndx_eval(irp_i,1,ispin) + i)
      call cont_eps_t3aaa(dcore,dcore(start_t3),
     &                    eps_ijk,irp_ijk,
     &                    ispin,
     &                    e4t,t3_ovrlp)
#endif /* _DO_EPS */

#ifdef _PROCESS_T3AAA
c   o pass the T3 amplitudes to a processing routine
      call process_t3aaa(
     &        dcore,dcore(start_scr),dim_scr,
     &        ispin,
     &        dcore(start_t3),max_abc,off_altbltc,
     &        i,irp_i,max_i,irp_ij,add_t2ij,
     &        j,irp_j,max_j,irp_jk,add_t2jk,
     &        k,irp_k,max_k,irp_ik,add_t2ik,irp_ijk,
     &        e5st
     &     )
#endif /* _PROCESS_T3AAA */

#ifndef ONLY_ONE_T3_DIST
c   o point to the next T3(IJK,*) distribution
      start_t3 = start_t3 + max_abc
#endif

c ----------------------------------------------------------------------

c     end if (mod(mpi_count,CNTOT).eq.CNID)
      end if

c      o increment and close the loop over orbitals I<J<K
            add_t2ij = add_t2ij + irpdpd(irp_ij,dist_altb)
            add_t2ik = add_t2ik + irpdpd(irp_ik,dist_altb)
            add_wij  = add_wij  + irpdpd(irp_ij,dist_ia)
            add_wik  = add_wik  + irpdpd(irp_ik,dist_ia)
c        end do i = 0, to_i
         end do
            add_t2jk = add_t2jk + irpdpd(irp_jk,dist_altb)
            add_wjk  = add_wjk  + irpdpd(irp_jk,dist_ia)
c        end do j = from_j, to_j
         end do
            if (irp_i.eq.irp_j) then
c            o if j starts at 1, instead of 0, skip another distribution
               add_t2jk = add_t2jk + irpdpd(irp_jk,dist_altb)
               add_wjk  = add_wjk  + irpdpd(irp_jk,dist_ia)
            end if
            if (irp_i.eq.irp_k) then
               start_t2ik =   start_t2ik
     &                      + ( irpdpd(irp_ik,dist_altb) * k )
               start_wik  =   start_wik
     &                      + ( irpdpd(irp_ik,dist_ia)   * k )
            else
               start_t2ik =   start_t2ik
     &                      + ( irpdpd(irp_ik,dist_altb) * max_i )
               start_wik  =   start_wik
     &                      + ( irpdpd(irp_ik,dist_ia)   * max_i )
            end if
c        end do k = from_k, max_k-1
         end do

c   o close the loop over irreps of occupieds in T3(IJK,ABC)
c     end if (do_irp_j)
      end if
c     end do irp_j = 1, irp_k
      end do
c     end if (max_k.ne.0)
      end if
c     end do irp_k = 1, nirrep
      end do

c     end if (max_abc.ne.0)
      end if
c     end do irp_ijk = 1, nirrep
      end do

      if (ispin.eq.1) then
         print *, 'There are ',mpi_count,' AAA triples.'
      else
         print *, 'There are ',mpi_count,' BBB triples.'
      end if

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

