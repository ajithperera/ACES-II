
c This routine contracts <cB|eI> with T2(kJ,eA) for all A, B, and c
c under a given I, J, and k.

c NOTE:
c    Orbitals I and J may be arbitrarily ordered. It is up to the
c calling routine to scale the contributions by -1 (via t3_factor)
c if J<I.

c INPUT
c    [yeah, right]

c ----------------------------------------------------------------------

#include "listnum.h"

c#define _DEBUG_CONT_GS_A2B1C2I1_T2AB

#ifdef _UNICOS
#  define B_GEMM sgemm
#else
#  define B_GEMM dgemm
#endif /* _UNICOS */

      subroutine cont_gs_a2b1c2i1_t2ab(
     &              dcore,  dscr,   dim_dscr,
     &              ispin1, ispin2,
     &              list_a2b1c2i1, dist_a1b2,
     &              i, irp_i, max_i,
     &              add_t2kj, off_a1b2,
     &              j, irp_j, max_j, irp_jk,
     &              k, irp_k, max_k, irp_ijk,
     &              t3_factor, start_t3, ndx_t3vrt
     &           )
      implicit none

c ARGUMENT LIST
      integer dim_dscr
      double precision dcore(*), dscr(dim_dscr)
      integer ispin1, ispin2
      integer list_a2b1c2i1, dist_a1b2
      integer i, irp_i, max_i
      integer add_t2kj, off_a1b2
      integer j, irp_j, max_j, irp_jk
      integer k, irp_k, max_k, irp_ijk
      double precision t3_factor
      integer start_t3, ndx_t3vrt(*)

c INTERNAL VARIABLES
      integer pq, irp_pq, max_pq
      integer e,  irp_e,  max_e, from_e,       last_e, num_e
      integer p,  irp_p,  max_p,         to_p
      integer q,  irp_q,  max_q
      integer r,  irp_r,  max_r, from_r, to_r, last_r, num_r
      logical do_irp_q, do_irp_r, do_t3
      integer size_scr, ioff, itmp
      integer ndx_ei, start_t2kjer, ndx, iadd, max_dists
      double precision dtmp_t3_factor

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"
#include "method.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      itmp = 0
c   o assert dim_dscr > 0
      if (dim_dscr.lt.1) then
         print *, '@CONT_GS_A2B1C2I1_T2AB: Assertion failed.'
         print *, '   dim_dscr = ',dim_dscr
         itmp = -1
      end if
c   o assert ispin1 and ispin2 are 1 or 2 and not equal
      if (.not.(((ispin1.eq.1).and.(ispin2.eq.2)).or.
     &          ((ispin1.eq.2).and.(ispin2.eq.1))    )) then
         print *, '@CONT_GS_A2B1C2I1_T2AB: Assertion failed.'
         print *, '   ispin1 = ',ispin1
         print *, '   ispin2 = ',ispin2
         itmp = -1
      end if
c   o assert i, j, and k are well-behaved
      if ((max_i.le.i).or.(i.lt.0).or.
     &    (max_j.le.j).or.(j.lt.0).or.
     &    (max_k.le.k).or.(k.lt.0)    ) then
         print *, '@CONT_GS_A2B1C2I1_T2AB: Assertion failed.'
         print *, '   i = ',i,', max_i = ',max_i
         print *, '   j = ',j,', max_j = ',max_j
         print *, '   k = ',k,', max_k = ',max_k
         itmp = -1
      end if
c   o assert irp_i, irp_j, and irp_k are in [1,nirrep]
      if ((irp_i  .lt.1).or.(nirrep.lt.irp_i  ).or.
     &    (irp_j  .lt.1).or.(nirrep.lt.irp_j  ).or.
     &    (irp_k  .lt.1).or.(nirrep.lt.irp_k  ).or.
     &    (irp_jk .lt.1).or.(nirrep.lt.irp_jk ).or.
     &    (irp_ijk.lt.1).or.(nirrep.lt.irp_ijk)    ) then
         print *, '@CONT_GS_A2B1C2I1_T2AB: Assertion failed.'
         print *, '   irp_i   = ',irp_i
         print *, '   irp_j   = ',irp_j
         print *, '   irp_k   = ',irp_k
         print *, '   irp_jk  = ',irp_jk
         print *, '   irp_ijk = ',irp_ijk
         itmp = -1
      end if
c   o assert dcore pointers are in [1,inf)
      if ((add_t2kj.lt.1).or.(start_t3.lt.1)) then
         print *, '@CONT_GS_A2B1C2I1_T2AB: Assertion failed.'
         print *, '   add_t2kj = ',add_t2kj
         print *, '   start_t3 = ',start_t3
         itmp = -1
      end if
      if (itmp.ne.0) then
         call xcc_errex
      end if
c   o warn if t3_factor = 0
      if (t3_factor.eq.0.0d0) then
         print *, '@CONT_GS_A2B1C2I1_T2AB: Warning.'
         print *, '   (t3_factor == 0.0d0) -> do nothing'
      end if
#endif /* _ASSERT */

      if (t3_factor.eq.0.0d0) return
c   o Even though add_t2kj points to the proper distribution,
c     load_t2_virts will negate the T2 amplitudes. We will reverse
c     this with t3_factor.
      if (ispin1.eq.1) then
         dtmp_t3_factor =  t3_factor
      else
         dtmp_t3_factor = -t3_factor
      end if

c The plan is:
c   o loop over irreps of R in T2(kJ,eR)
c      o loop over orbital e
c         - load as many <pQ|eI> lists as possible (all pQ for each e)
c         - if (ispin1.eq.1) then <Qp|Ie> -> <pQ|eI>
c         - <pQ|eI> -> <Pq|eI> (to use put_t3aba.F.src)
c         o loop over orbital R in T2(kJ,eR)
c            - create T2(kJ,eR) (all e for each R)
c            - gemm <Pq|eI> with T2(kJ,eR)
c            o loop over irreps of q in T3(IJk,PqR)
c               - add in the T3(IJk,PqR) contributions

c   o loop over irreps of R in T2(kJ,eR)
      do irp_r = 1, nirrep
c      o check for T2 amplitudes
         irp_e = dirprd(irp_r,irp_jk)
         max_r = vrt(irp_r,ispin1)
         max_e = vrt(irp_e,ispin2)
         do_irp_r = ((max_e.ne.0).and.(max_r.ne.0))
         if (do_irp_r) then
c         o check for integrals
            irp_pq = dirprd(irp_e,irp_i)
            max_pq = irpdpd(irp_pq,dist_a1b2)
            do_irp_r = (max_pq.ne.0)
         end if
      if (do_irp_r) then

c      o initialize values for these irreps
         last_e = max_e-1
         last_r = max_r-1
c        size_scr : scratch space minus at least one T3(IJk,PqR) product
c                   matrix and T2(kJ,eR) factor vector
         size_scr = dim_dscr - ( max_pq + max_e )
         itmp = mod(size_scr,max_pq)
         max_dists = ( size_scr - itmp ) / max_pq
         if (max_dists.eq.0) then
            print *, '@CONT_GS_A2B1C2I1_T2AB: Not enough memory.'
            print *, '                        At least ',
     &               lshift((max_pq-size_scr),3),
     &               ' more bytes are needed.'
            call xcc_errex
         end if

c      o point to the first |Ie> (or |Ei>) list
         if ((ispin1.eq.1).and.(iuhf.eq.1)) then
            itmp = nirrep*(irp_e-1) + irp_i
            ndx_ei = 1 + off_pq(itmp,h_Offsets_Ia) + i
         else
            itmp = nirrep*(irp_i-1) + irp_e
            ndx_ei = 1 + off_pq(itmp,h_Offsets_Ai) + ( max_e * i )
         end if

c      o loop over orbital e in <pQ|eI>
         from_e = 0
         do while (from_e.lt.max_e)

c         o get as many <Pq|eI> lists as possible
            num_e = max_e-from_e
            num_e = min(num_e,max_dists)
            if (iuhf.eq.0) then
c            o RHF can be "hardcoded"
               call getlst(dscr,ndx_ei,num_e,2,irp_pq,h_GS_AbCi_Ab_Ci)
               call transp_dist(dscr,vrt(1,1),vrt(1,1),
     &                          irp_pq,num_e,.false.)
               ndx_ei = ndx_ei + num_e
            else
               if (ispin1.eq.1) then
                  do e = 0, num_e-1
                     call getlst(dscr(1+(e*max_pq)),ndx_ei,1,2,irp_pq,
     &                           list_a2b1c2i1)
                     ndx_ei = ndx_ei + max_i
                  end do
               else
                  call getlst(dscr(1),ndx_ei,num_e,2,irp_pq,
     &                        list_a2b1c2i1)
                  call transp_dist(dscr(1),vrt(1,1),vrt(1,2),
     &                             irp_pq,num_e,.false.)
                  ndx_ei = ndx_ei + num_e
               end if
c           end if (iuhf.eq.0)
            end if

c         o reset the T2(kJ,eR) pointer
c           This is the scratch address of the temporary factor.
            start_t2kjer = 1 + ( num_e * max_pq )

c         o loop over orbital R in T2(kJ,eR)
            from_r = 0
            do while (from_r.lt.max_r)

c            o We know there is enough memory to process at least one R orbital
c              (from size_scr), so we do not need to check for memory.
c              ndx : (temp) a T2, T3 column pair (factor, product)
               ndx   = max_e + max_pq
               itmp  = 1 + dim_dscr - start_t2kjer
               ioff  = mod(itmp,ndx)
c              to_r : (temp) the number of column pairs that can fit in dscr
               to_r  = ( itmp - ioff ) / ndx
               num_r = min(to_r,(max_r-from_r))
               to_r  = from_r + num_r - 1

c            o create T2(kJ,eR) (all e for each R)
               call load_t2_virts(
     &                 dcore(add_t2kj), dscr(start_t2kjer),
     &                 irp_e, max_e,               ispin2,
     &                 irp_r, max_r, from_r, to_r, ispin1, off_a1b2
     &              )

c            o gemm <Pq|eI> with T2(kJ,eR)
               iadd = start_t2kjer + from_e
               ndx  = start_t2kjer + ( max_e * num_r )
               call B_GEMM('N','N',max_pq,num_r,num_e,
     &                     dtmp_t3_factor,dscr(1),   max_pq,
     &                                    dscr(iadd),max_e,
     &                     0.0d0,         dscr(ndx), max_pq)
c      if (ispin1.eq.1) then
c      print *,
c      print *, '@CONT_GS_A2B1C2I1_T2AB: <Pq|eI>: (',num_e,' cols)'
c      call dmat_prt(dscr(1),max_pq,num_e,max_pq)
c      print *, '@CONT_GS_A2B1C2I1_T2AB: T2(kJ,eR): (',num_r,' cols)'
c      call dmat_prt(dscr(iadd),num_e,num_r,max_e)
c      print *,
c      end if


c            o loop over T3(IJk,PqR) contributions
               do r = from_r, to_r

c ----------------------------------------------------------------------

c   o add in the T3 contributions
#include "put_t3aba.F.src"

c ----------------------------------------------------------------------

c            o close the loop over T3(IJk,PqR) contributions
c              end do r = from_r, to_r
               end do

c         o increment and close the loop over orbital R in T2(kJ,eR)
               from_r = from_r + num_r
c           end do while (from_r.lt.max_r)
            end do

c      o increment and close the loop over orbital e in <Pq|eI>
            from_e = from_e + num_e
c        end do while (from_e.lt.max_e)
         end do

c   o close the loop over irreps of R in T2(kJ,eR)
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

