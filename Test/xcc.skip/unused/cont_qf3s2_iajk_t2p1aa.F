
c This routine contracts T2(1)(MI,AB) with <MC{}JK> for all A, B, and C
c under a given I, J, and K.

c NOTE:
c    Orbitals I, J, and K may be arbitrarily ordered. It is up to the
c calling routine to scale the contributions by -1 (via s3_factor)
c if J<I<K.

c INPUT
c    [yeah, right]

c ----------------------------------------------------------------------

c#define _DEBUG_CONT_QF3S2_IAJK_T2P1AA

      subroutine cont_qf3s2_iajk_t2p1aa(
     &              dcore, dscr, dim_dscr,
     &              ispin,
     &              off_iltj, dist_altb, dist_ia,
     &              i, irp_i, max_i, irp_ij,
     &              j, irp_j, max_j, irp_jk, ioff_jk,
     &              k, irp_k, max_k, irp_ijk,
     &              s3_factor, s3, off_altbltc
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr), s3(*)
      integer ispin, off_iltj, dist_altb, dist_ia, off_altbltc(*)
      integer i, irp_i, max_i, irp_ij
      integer j, irp_j, max_j, irp_jk, ioff_jk
      integer k, irp_k, max_k, irp_ijk
      double precision s3_factor

c INTERNAL VARIABLES
      integer irp_pq, max_pq
      integer irp_m,  max_m
      integer irp_r,  max_r
      logical do_irp_r
      integer itmp
      integer start_scr, size_scr
      integer start_s2, add_s2
      integer start_t2p1

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      itmp = 0
c   o assert dim_dscr > 0
      if (dim_dscr.lt.1) then
         write(*,*) '@CONT_QF3S2_IAJK_T2P1AA: Assertion failed.'
         write(*,*) '   dim_dscr = ',dim_dscr
         itmp = -1
      end if
c   o assert ispin is 1 or 2
      if ((ispin.lt.1).or.(2.lt.ispin)) then
         write(*,*) '@CONT_QF3S2_IAJK_T2P1AA: Assertion failed.'
         write(*,*) '   ispin = ',ispin
         itmp = -1
      end if
c   o assert off_iltj, dist_altb, and dist_ia are all natural
      if ((off_iltj.lt.1).or.(dist_altb.lt.1).or.(dist_ia.lt.1)) then
         write(*,*) '@CONT_QF3S2_IAJK_T2P1AA: Assertion failed.'
         write(*,*) '   off_iltj  = ',off_iltj
         write(*,*) '   dist_ia   = ',dist_ia
         write(*,*) '   dist_altb = ',dist_altb
         itmp = -1
      end if
c   o assert i, j, and k are well-behaved
      if ((max_i.le.i).or.(i.lt.0).or.
     &    (max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         write(*,*) '@CONT_QF3S2_IAJK_T2P1AA: Assertion failed.'
         write(*,*) '   i = ',i,', max_i = ',max_i
         write(*,*) '   j = ',j,', max_j = ',max_j
         write(*,*) '   k = ',k,', max_k = ',max_k
         itmp = -1
      end if
c   o assert irp_i, irp_j, and irp_k are in [1,nirrep]
      if ((irp_i  .lt.1).or.(nirrep.lt.irp_i  ).or.
     &    (irp_j  .lt.1).or.(nirrep.lt.irp_j  ).or.
     &    (irp_k  .lt.1).or.(nirrep.lt.irp_k  ).or.
     &    (irp_ij .lt.1).or.(nirrep.lt.irp_ij ).or.
     &    (irp_jk .lt.1).or.(nirrep.lt.irp_jk ).or.
     &    (irp_ijk.lt.1).or.(nirrep.lt.irp_ijk)    ) then
         write(*,*) '@CONT_QF3S2_IAJK_T2P1AA: Assertion failed.'
         write(*,*) '   irp_i   = ',irp_i
         write(*,*) '   irp_j   = ',irp_j
         write(*,*) '   irp_k   = ',irp_k
         write(*,*) '   irp_ij  = ',irp_ij
         write(*,*) '   irp_jk  = ',irp_jk
         write(*,*) '   irp_ijk = ',irp_ijk
         itmp = -1
      end if
c   o assert ioff_jk >= 0
      if (ioff_jk.lt.0) then
         write(*,*) '@CONT_QF3S2_IAJK_T2P1AA: Assertion failed.'
         write(*,*) '   ioff_jk = ',ioff_jk
         itmp = -1
      end if
      if (itmp.ne.0) then
         call xcc_errex
      end if
c   o warn if s3_factor = 0
      if (s3_factor.eq.0.0d0) then
         write(*,*) '@CONT_QF3S2_IAJK_T2P1AA: Warning.'
         write(*,*) '   (s3_factor == 0.0d0) -> do nothing'
      end if
#endif /* _ASSERT */

      if (s3_factor.eq.0.0d0) return

c ----------------------------------------------------------------------

c The plan is:
c   - get the <MR{}JK> distribution
c   o loop over irreps of R in <MR{}JK>
c      - load/unpack T2(1)(MI,PQ) (all P<Q for all M)
c      o loop over batches of R orbitals in <MR{}JK>
c         - contract T2(1)(MI,PQ) with <MR{}JK> and add into S3

c   o initialize values
      start_scr = 1
      size_scr  = dim_dscr

c   o create the <MR{}JK> distribution
      max_pq = irpdpd(irp_jk,dist_ia)
      if (max_pq.eq.0) return
      start_s2  = start_scr
      start_scr = start_scr + max_pq
      size_scr  = size_scr  - max_pq
      if (size_scr.lt.0) then
         write(*,*) '@CONT_QF3S2_IAJK_T2P1AA: Not enough memory.'
         write(*,*) '                         At least ',
     &              lshift(-size_scr,3),
     &              ' more bytes are needed.'
         call xcc_errex
      end if
      call make_qf3_hh_aijk(
     &        dcore, dscr(start_scr), size_scr,
     &        ispin, dscr(start_s2),
     &        j, irp_j, max_j, irp_jk, ioff_jk,
     &        k, irp_k, max_k, irp_jk
     &     )
c   o THESE ARE NOT NEGATED; therefore, negate gemm.
      call transp_dist(
     &        dscr(start_s2),
     &        vrt(1,ispin), pop(1,ispin),
     &        irp_jk, 1,
     &        .false.
     &     )

c   o loop over irreps of R in <MR{}JK>
      start_t2p1 = start_scr
      add_s2 = start_s2
      do irp_r = 1, nirrep
c      o check for integrals
         irp_m = dirprd(irp_r,irp_jk)
         max_m = pop(irp_m,ispin)
         max_r = vrt(irp_r,ispin)
         do_irp_r = ((max_r.ne.0).and.(max_m.ne.0))
         if (do_irp_r) then
c         o check for T2(1) amplitudes
            irp_pq = dirprd(irp_m,irp_i)
            max_pq = irpdpd(irp_pq,dist_altb)
            do_irp_r = ((max_pq.ne.0).and.
     &                  (((irp_pq.eq.1).and.(max_m.ne.1)).or.
     &                   (irp_pq.ne.1)                       )
     &                 )
         end if
      if (do_irp_r) then

c      o Make sure there is enough room for all T2(1)(MI,PQ) and at least
c        one S3(IJK,PQR) product. This should be changed in the future
c        to allow for batch processing of M orbitals.
         start_scr = start_t2p1 + ( max_pq * max_m )
         size_scr  = 1 + dim_dscr - start_scr
         itmp = size_scr - max_pq
         if (itmp.lt.0) then
            write(*,*) '@CONT_QF3S2_IAJK_T2P1AA: Not enough memory.'
            write(*,*) '                         At least ',
     &                 lshift(-itmp,3),
     &                 ' more bytes are needed.'
            call xcc_errex
         end if

c      o load/unpack T2(1)(MI,PQ) (all P<Q for all M)
         if (irp_m.lt.irp_i) then
            itmp = rshift(irp_i*(irp_i-1),1) + irp_m
         else
            itmp = rshift(irp_m*(irp_m-1),1) + irp_i
         end if
         itmp = off_pltp(itmp,off_iltj)
         itmp = ndx_t2p1(irp_pq,ispin) + ( max_pq * itmp )
         call extract_s2_dists(
     &           dcore(itmp),  dscr(start_t2p1),  max_pq,
     &           irp_m, max_m, 0, max_m-1, max_m,
     &           irp_i, max_i, i,
     &           .true., .false.
     &        )

c ----------------------------------------------------------------------

c   o s3_factor is negative because the intermediates were not negated.
      call cont_s3pltpltp_s2pltp_s2xp(
     &        dscr(start_t2p1), irp_pq, max_pq, max_m,
     &        dscr(add_s2),     irp_r,  max_m,  0,     max_r-1, max_r,
     &        -s3_factor, s3, off_altbltc,
     &        vrt(1,ispin), dscr(start_scr), size_scr
     &     )

c ----------------------------------------------------------------------

c   o close the loop over irreps of R in <MR{}JK>
c     end if (do_irp_r)
      end if
         add_s2 = add_s2 + ( max_m * max_r )
c     end do irp_r = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

