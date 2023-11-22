
c This routine contracts T2(Mk,Bc) with <MA||IJ> for all A, B, and c
c under a given I, J, and k.

c NOTE:
c    Orbital I must be less than orbital J for the antisymmetry to
c wash out; therefore, unless you know what you are doing, t3_factor
c should be positive.

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

c#define _DEBUG_CONT_HH_IAJK_T2AB

#ifdef _UNICOS
#  define B_GEMM sgemm
#else
#  define B_GEMM dgemm
#endif /* _UNICOS */

      subroutine cont_hh_iajk_t2ab(
     &              dcore,  dscr,   dim_dscr,
     &              ispin1, ispin2,
     &              off_i1j2, dist_a1b2,
     &              k, irp_k, max_k,
     &              add_wij,
     &              i, irp_i, max_i, irp_ij,
     &              j, irp_j, max_j, irp_ijk,
     &              t3_factor, start_t3, ndx_t3vrt
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr)
      integer ispin1, ispin2
      integer off_i1j2, dist_a1b2
      integer k, irp_k, max_k
      integer add_wij
      integer i, irp_i, max_i, irp_ij
      integer j, irp_j, max_j, irp_ijk
      double precision t3_factor
      integer start_t3, ndx_t3vrt(*)

c INTERNAL VARIABLES
      integer pq, irp_pq, max_pq
      integer     irp_m,  max_m
      integer p,  irp_p,  max_p,         to_p
      integer q,  irp_q,  max_q
      integer r,  irp_r,  max_r, from_r, to_r, num_r
      logical do_irp_r,       do_irp_q,       do_t3
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
         print *, '@CONT_HH_IAJK_T2AB: Assertion failed.'
         print *, '   dim_dscr = ',dim_dscr
         itmp = 1
      end if
c   o assert ispin1 and ispin2 are 1 or 2 and not equal
      if (.not.(((ispin1.eq.1).and.(ispin2.eq.2)).or.
     &          ((ispin1.eq.2).and.(ispin2.eq.1))    )) then
         print *, '@CONT_HH_IAJK_T2AB: Assertion failed.'
         print *, '   ispin1 = ',ispin1
         print *, '   ispin2 = ',ispin2
         itmp = 1
      end if
c   o assert i, j, and k are well-behaved
      if ((max_i.le.i).or.(i.lt.0).or.
     &    (max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         print *, '@CONT_HH_IAJK_T2AB: Assertion failed.'
         print *, '   i = ',i,', max_i = ',max_i
         print *, '   j = ',j,', max_j = ',max_j
         print *, '   k = ',k,', max_k = ',max_k
         itmp = 1
      end if
c   o assert i<j
      if ((irp_i.gt.irp_j).or.((irp_i.eq.irp_j).and.(i.gt.j))) then
         print *, '@CONT_HH_IAJK_T2AB: Assertion failed.'
         print *, '   i = ',1+i,'[',irp_i,']'
         print *, '   j = ',1+j,'[',irp_j,']'
         itmp = 1
      end if
c   o assert irp_i, irp_j, and irp_k are in [1,nirrep]
      if ((irp_i  .lt.1).or.(nirrep.lt.irp_i  ).or.
     &    (irp_j  .lt.1).or.(nirrep.lt.irp_j  ).or.
     &    (irp_k  .lt.1).or.(nirrep.lt.irp_k  ).or.
     &    (irp_ij .lt.1).or.(nirrep.lt.irp_ij ).or.
     &    (irp_ijk.lt.1).or.(nirrep.lt.irp_ijk)    ) then
         print *, '@CONT_HH_IAJK_T2AB: Assertion failed.'
         print *, '   irp_i   = ',irp_i
         print *, '   irp_j   = ',irp_j
         print *, '   irp_k   = ',irp_k
         print *, '   irp_ij  = ',irp_ij
         print *, '   irp_ijk = ',irp_ijk
         itmp = 1
      end if
c   o assert dcore pointers are in [1,inf)
      if ((add_wij.lt.1).or.(start_t3.lt.1)) then
         print *, '@CONT_HH_IAJK_T2AB: Assertion failed.'
         print *, '   add_wij  = ',add_wij
         print *, '   start_t3 = ',start_t3
         itmp = 1
      end if
      if (itmp.ne.0) call xcc_errex
c   o warn if t3_factor = 0
      if (t3_factor.eq.0.0d0) then
         print *, '@CONT_HH_IAJK_T2AB: Warning.'
         print *, '   (t3_factor == 0.0d0) -> do nothing'
      end if
#endif /* _ASSERT */

      if (t3_factor.eq.0.0d0) return

c The plan is:
c   o loop over irreps of R in <MR||IJ>
c      - create T2(Mk,Pq) (all Pq for all M)
c      o loop over orbital R in <MR||IJ>
c         - gemm T2(Mk,Pq) with <MR||IJ>
c         o loop over irreps of q in T3(IJk,PqR)
c            - add in the T3(IJk,PqR) contributions

      start_w = add_wij

c   o loop over irreps of R in <MR||IJ>
      do irp_r = 1, nirrep
c      o check for integrals
         irp_m = dirprd(irp_r,irp_ij)
         max_m = pop(irp_m,ispin1)
         max_r = vrt(irp_r,ispin1)
         do_irp_r = ((max_r.ne.0).and.(max_m.ne.0))
         if (do_irp_r) then
c         o check for T2 amplitudes
            irp_pq = dirprd(irp_m,irp_k)
            max_pq = irpdpd(irp_pq,dist_a1b2)
            do_irp_r = (max_pq.ne.0)
         end if
      if (do_irp_r) then

c      o Make sure there is enough room for all T2(Mk,Pq) and at least
c        one T3(IJk,PqR) product. This should be changed in the future
c        to allow for batch processing of M orbitals.
         start_product = 1 + ( max_pq * max_m )
         size_scr = 1 + dim_dscr - start_product
         itmp = mod(size_scr,max_pq)
         max_dists = ( size_scr - itmp ) / max_pq
         if (max_dists.eq.0) then
            print *, '@CONT_HH_IAJK_T2AB: Not enough memory.'
            print *, '                    At least ',
     &               lshift((max_pq-size_scr),3),
     &               ' more bytes are needed.'
            call xcc_errex
         end if

c      o create T2(Mk,Pq) (all Pq for all M)
         call load_t2_dists(
     &           dcore, dscr,
     &           irp_pq, max_pq,
     &           irp_m,  max_m,  0, max_m-1, ispin1,
     &           irp_k,  max_k,  k,          ispin2, off_i1j2
     &        )

c      o loop over orbital R in <MR||IJ>
         from_r = 0
         do while (from_r.lt.max_r)

            num_r = max_r-from_r
            num_r = min(num_r,max_dists)
            to_r  = from_r + num_r - 1

            ndx = start_product

c         o gemm T2(Mk,Pq) with <MR||IJ>
            call B_GEMM('N','N',max_pq,num_r,max_m,
     &                  t3_factor,dscr(1),       max_pq,
     &                            dcore(start_w),max_m,
     &                  0.0d0,    dscr(ndx),     max_pq)
c      print *,
c      print *, '@D1T3: T2(MK,PQ): (',max_m,' cols)'
c      call dmat_prt(dscr,max_pq,max_m,max_pq)
c      print *, '@D1T3: <MR||IJ>: (',num_r,' cols)'
c      call dmat_prt(dcore(start_w),max_m,num_r,max_m)
c      print *,

c         o point to the next batch of <MR||IJ> integrals
            start_w = start_w + ( max_m * num_r )

c         o loop over orbital R in T3(IJk,PqR)
            do r = from_r, to_r

c ----------------------------------------------------------------------

c   o add in the T3 contributions
#include "put_t3aba.F.src"

c ----------------------------------------------------------------------

c         o close the loop over orbital R in T3(IJk,PqR)
c           end do r = from_r, to_r
            end do

c      o increment and close the loop over orbital R in <MR||IJ>
            from_r = from_r + num_r
c        end do while (from_r.lt.max_r)
         end do

c   o close the loop over irreps of R in <MR||IJ>
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

