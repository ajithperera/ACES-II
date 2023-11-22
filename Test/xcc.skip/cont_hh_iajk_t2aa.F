
c This routine contracts T2(MI,AB) with <MC||JK> for all A, B, and C
c under a given I, J, and K.

c NOTE:
c    Orbitals I, J, and K may be arbitrarily ordered. It is up to the
c calling routine to scale the contributions by -1 (via t3_factor)
c if J<I<K.

c NOTE:
c    If you are unlucky enough to have to debug this, then there are a
c few things which you need to know...
c 1) Pay CLOSE attention to the scope of each variable. itmp, iadd,
c    ndx, ioff, etc. are used all over the place as temporary
c    variables. If this were C (instead of Fortran), then each segment
c    needing a temp integer would be surrounded by braces.
c 2) [more to come as the code gets more complex]

c INPUT
c    [yeah, right]

c ----------------------------------------------------------------------

c#define _DEBUG_CONT_HH_IAJK_T2AA

#ifdef _UNICOS
#  define B_GEMM sgemm
#else
#  define B_GEMM dgemm
#endif /* _UNICOS */

      subroutine cont_hh_iajk_t2aa(
     &              dcore, dscr, dim_dscr,
     &              ispin,
     &              off_iltj, dist_altb,
     &              i, irp_i, max_i,
     &              add_wjk,
     &              j, irp_j, max_j, irp_jk,
     &              k, irp_k, max_k, irp_ijk,
     &              t3_factor, start_t3, ndx_t3vrt
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr)
      integer ispin, off_iltj, dist_altb, add_wjk
      integer i, irp_i, max_i
      integer j, irp_j, max_j, irp_jk
      integer k, irp_k, max_k, irp_ijk
      double precision t3_factor
      integer start_t3, ndx_t3vrt(*)

c INTERNAL VARIABLES
      integer pq, irp_pq, max_pq
      integer     irp_m,  max_m
      integer p,  irp_p,  max_p,         to_p
      integer q,  irp_q,  max_q,         to_q
      integer r,  irp_r,  max_r, from_r, to_r, num_r
      logical do_irp_r, do_irp_q, do_t3
      integer size_scr, ioff, itmp
      integer start_w, start_product, iadd, ndx, max_dists

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      itmp = 0
c   o assert dim_dscr > 0
      if (dim_dscr.lt.1) then
         print *, '@CONT_HH_IAJK_T2AA: Assertion failed.'
         print *, '   dim_dscr = ',dim_dscr
         itmp = 1
      end if
c   o assert ispin is 1 or 2
      if ((ispin.lt.1).or.(2.lt.ispin)) then
         print *, '@CONT_HH_IAJK_T2AA: Assertion failed.'
         print *, '   ispin = ',ispin
         itmp = 1
      end if
c   o assert something about dist_altb?
c   o assert i, j, and k are well-behaved
      if ((max_i.le.i).or.(i.lt.0).or.
     &    (max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         print *, '@CONT_HH_IAJK_T2AA: Assertion failed.'
         print *, '   i = ',i,', max_i = ',max_i
         print *, '   j = ',j,', max_j = ',max_j
         print *, '   k = ',k,', max_k = ',max_k
         itmp = 1
      end if
c   o assert irp_i, irp_j, and irp_k are in [1,nirrep]
      if ((irp_i  .lt.1).or.(nirrep.lt.irp_i  ).or.
     &    (irp_j  .lt.1).or.(nirrep.lt.irp_j  ).or.
     &    (irp_k  .lt.1).or.(nirrep.lt.irp_k  ).or.
     &    (irp_jk .lt.1).or.(nirrep.lt.irp_jk ).or.
     &    (irp_ijk.lt.1).or.(nirrep.lt.irp_ijk)    ) then
         print *, '@CONT_HH_IAJK_T2AA: Assertion failed.'
         print *, '   irp_i   = ',irp_i
         print *, '   irp_j   = ',irp_j
         print *, '   irp_k   = ',irp_k
         print *, '   irp_jk  = ',irp_jk
         print *, '   irp_ijk = ',irp_ijk
         itmp = 1
      end if
c   o assert dcore pointers are in [1,inf)
      if ((add_wjk.lt.1).or.(start_t3.lt.1)) then
         print *, '@CONT_HH_IAJK_T2AA: Assertion failed.'
         print *, '   add_wjk  = ',add_wjk
         print *, '   start_t3 = ',start_t3
         itmp = 1
      end if
      if (itmp.ne.0) call xcc_errex
c   o warn if t3_factor = 0
      if (t3_factor.eq.0.0d0) then
         print *, '@CONT_HH_IAJK_T2AA: Warning.'
         print *, '   (t3_factor == 0.0d0) -> do nothing'
      end if
#endif /* _ASSERT */

      if (t3_factor.eq.0.0d0) return

c The plan is:
c   o loop over irreps of R in <MR||JK>
c      - create T2(MI,PQ) (all P<Q for all M)
c      o loop over orbital R in <MR||JK>
c         - gemm T2(MI,PQ) with <MR||JK>
c         o loop over irreps of Q in T3(IJK,PQR)
c            - add in the T3(IJK,PQR) contributions

      start_w = add_wjk

c   o loop over irreps of R in <MR||JK>
      do irp_r = 1, nirrep
c      o check for integrals
         irp_m = dirprd(irp_r,irp_jk)
         max_m = pop(irp_m,ispin)
         max_r = vrt(irp_r,ispin)
         do_irp_r = ((max_r.ne.0).and.(max_m.ne.0))
         if (do_irp_r) then
c         o check for T2 amplitudes
            irp_pq = dirprd(irp_m,irp_i)
            max_pq = irpdpd(irp_pq,dist_altb)
            do_irp_r = ((max_pq.ne.0).and.
     &                  (((irp_pq.eq.1).and.(max_m.ne.1)).or.
     &                   (irp_pq.ne.1)                       )
     &                 )
         end if
      if (do_irp_r) then

c      o Make sure there is enough room for all T2(MI,PQ) and at least
c        one T3(IJK,PQR) product. This should be changed in the future
c        to allow for batch processing of M orbitals.
         start_product = 1 + ( max_pq * max_m )
         size_scr = 1 + dim_dscr - start_product
         itmp = mod(size_scr,max_pq)
         max_dists = ( size_scr - itmp ) / max_pq
         if (max_dists.eq.0) then
            print *, '@CONT_HH_IAJK_T2AA: Not enough memory.'
            print *, '                    At least ',
     &               lshift((max_pq-size_scr),3),
     &               ' more bytes are needed.'
            call xcc_errex
         end if

c      o create T2(MI,PQ) (all P<Q for all M)
         call load_t2_dists(
     &           dcore, dscr,
     &           irp_pq, max_pq,
     &           irp_m,  max_m,  0, max_m-1, ispin,
     &           irp_i,  max_i,  i,          ispin, off_iltj
     &        )

c      o loop over orbital R in <MR||JK>
         from_r = 0
         do while (from_r.lt.max_r)

            num_r = max_r-from_r
            num_r = min(num_r,max_dists)
            to_r  = from_r + num_r - 1

            ndx = start_product

c         o gemm T2(MI,PQ) with <MR||JK>
            call B_GEMM('N','N',max_pq,num_r,max_m,
     &                  t3_factor,dscr(1),       max_pq,
     &                            dcore(start_w),max_m,
     &                  0.0d0,    dscr(ndx),     max_pq)
c      print *,
c      print *, '@CONT_HH_IAJK_T2AA: T2(MI,PQ): (',max_m,' cols)'
c      call dmat_prt(dscr,max_pq,max_m,max_pq)
c      print *, '@CONT_HH_IAJK_T2AA: <MR||JK>: (',num_r,' cols)'
c      call dmat_prt(dcore(start_w),max_m,num_r,max_m)
c      print *,

c         o point to the next batch of <MR||JK> integrals
            start_w = start_w + ( max_m * num_r )

c         o loop over orbital R in T3(IJK,PQR)
            do r = from_r, to_r

c ----------------------------------------------------------------------

c   o add in the T3 contributions
#include "put_t3aaa.F.src"

c ----------------------------------------------------------------------

c         o close the loop over orbital R in T3(IJK,PQR)
c           end do r = from_r, to_r
            end do

c      o increment and close the loop over orbital R in <MR||JK>
            from_r = from_r + num_r
c        end do while (from_r.lt.max_r)
         end do

c   o close the loop over irreps of R in <MR||JK>
c     else if (.not.do_irp_r) then
      else
c      o increment start_w even if there are no T2 amplitudes
         start_w = start_w + ( max_m * max_r )
c     end if (do_irp_r)
      end if
c     end do irp_r = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

