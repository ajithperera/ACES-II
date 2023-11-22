
c This routine creates the antisymmetric product of a 1-particle operator
c squared. The 2-particle result is scaled and added into the packed 2-particle
c destination. The symmetry-forbidden elements of the product are not touched.

c The core of the subroutine performs the following:
c    S(p<q,r<s) += alpha * [ S(p,r)*S(q,s) - S(q,r)*S(p,s) ]

c Since S1 is squared, S2 must always be totally symmetric. This routine further
c assumes the destination is packed.

c INPUT
c int irp_rs   : the total ket symmetry of the S2 operator
c int rs_start : the rs super-index to start on (think getlst columns)
c int num_rs   : the number of rs pairs to update
c double alpha : the S1*S1 scaling factor
c double s1(*) : the S1 matrix elements for all ket irreps
c int irp_s1   : the total symmetry of the S1 operator
c int s1_bra_pop(*) : the population vector for the S1 bra elements
c int s1_ket_pop(*) : the population vector for the S1 ket elements
c int off_s1(*) : the vector of offsets to reach each S1 ket irrep

c INPUT/OUTPUT
c double s2(*) : the array to update with S1*S1

c#define _DEBUG_ANTISYM_S2PACK_S1PQ_SQR
c#define _USE_XGEMM

      subroutine antisym_s2pack_s1pq_sqr(
     &                    s2,irp_rs,rs_start,num_rs,
     &              alpha,s1,irp_s1,s1_bra_pop,s1_ket_pop,off_s1
     &           )
      implicit none

c ARGUMENT LIST
      double precision alpha, s1(*), s2(*)
      integer irp_rs, rs_start, num_rs
      integer irp_s1, s1_bra_pop(*), s1_ket_pop(*), off_s1(*)

c EXTERNAL FUNCTIONS
      integer rshift

c INTERNAL VARIABLES
      double precision dTmp_r, dTmp_s
      integer s1r_ptr2, s1r_ptr, s1s_ptr, s2_ndx
      integer rs_lo, rs, rs_hi
      integer irp_s1q, max_s1q
      integer irp_r, r, max_r, to_r
      integer irp_s, s, max_s
      integer irp_p, p, max_p, to_p
      integer irp_q, q, max_q
      logical do_irp_s, do_irp_q
      integer iTmp

c COMMON BLOCKS
#include "syminf.com" /* for nirrep and dirprd(,) */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert irp_rs and irp_s1 are in [1,nirrep]
      if ((irp_rs.lt.1).or.(nirrep.lt.irp_rs).or.
     &    (irp_s1.lt.1).or.(nirrep.lt.irp_s1)    ) then
         print *, '@ANTISYM_S2PACK_S1PQ_SQR: Assertion failed.'
         print *, '   irp_rs = ',irp_rs
         print *, '   irp_s1 = ',irp_s1
         iTmp = 1
      end if
c   o assert rs_start and num_rs are properly bound
      if ((rs_start.lt.1).or.(num_rs.lt.0)) then
         print *, '@ANTISYM_S2PACK_S1PQ_SQR: Assertion failed.'
         print *, '   rs_start = ',rs_start
         print *, '   num_rs   = ',num_rs
         iTmp = 1
      end if
      if (iTmp.ne.0) call aces_exit(iTmp)
#endif

      if ((num_rs.lt.1).or.(alpha.eq.0.d0)) return

c ----------------------------------------------------------------------

c   o initialize the various counters, indices, boundaries, and pointers
      s2_ndx = 1
      rs_lo  = rs_start - 1
      rs_hi  = rs_start + num_rs
      rs     = 1

c NOTE:
c    We start counting the S2 columns from the beginning regardless of
c whether rs_start is 1 or not. The reason is that the alternative of
c starting exactly at rs_start would add code equally, if not more,
c complex than what we have now. Additionally, the amount of work
c performed inside an rs-pair loop dwarfs the expense of the conditional.

c ----------------------------------------------------------------------

      if (irp_rs.eq.1) then

c   o loop over the S2 ket irreps
      do irp_s = 1, nirrep
         max_s = s1_ket_pop(irp_s)
      if (max_s.gt.1) then

c      o for these ket irreps
c         - cache the direct product look-up for the symmetry-allowed bra irrep
c         - cache the population for the bra irrep
c         - start the S1 column pointers
         irp_s1q = dirprd(irp_s,irp_s1)
         max_s1q = s1_bra_pop(irp_s1q)
         s1r_ptr2 = 1 + off_s1(irp_s)
         s1s_ptr  = s1r_ptr2 + max_s1q

c      o loop over the r<s orbital pairs (as offsets)
         do s = 1, max_s-1
c         o reset the S1(P,R) pointer
            s1r_ptr = s1r_ptr2
         do r = 0, s-1

c         o do work if rs is above the low bound
            if (rs.gt.rs_lo) then

c            o loop over the S2 bra irreps
               do irp_q = 1, nirrep
                  max_q = s1_bra_pop(irp_q)
               if (max_q.gt.1) then

c ----------------------------------------------------------------------

c   o process symmetry-allowed elements
      if (irp_q.eq.irp_s1q) then

c      o S2(pq,rs) += alpha*S1(p,r)*S1(q,s) - alpha*S1(q,r)*S1(p,s)
         do q = 1, max_q-1
            dTmp_s =  alpha*s1(s1s_ptr+q)
            dTmp_r = -alpha*s1(s1r_ptr+q)
            do p = 0, q-1
               s2(s2_ndx+p) = s2(s2_ndx+p) + dTmp_s*s1(s1r_ptr+p)
     &                                     + dTmp_r*s1(s1s_ptr+p)
            end do
            s2_ndx = s2_ndx + q
         end do

c     else if (Q != S x S1) then
      else

c      o point to the next set of S2 p<q pairs
         s2_ndx = s2_ndx + rshift(max_q*(max_q-1),1)

c   o close the symmetry conditional
c     end if ((P == R x S1) && (Q == S x S1))
      end if

c ----------------------------------------------------------------------

c            o close the loop over the S2 bra irreps
c              end if (max_q.gt.1)
               end if
c              end do irp_q = 1, nirrep
               end do

c         o close the low bound conditional
c           end if (rs.gt.rs_lo)
            end if

c         o increment the rs super-index and quit if necessary
            rs = rs + 1
            if (rs.eq.rs_hi) return

c      o increment the S1 pointers and close the loop over the r<s orbital pairs
            s1r_ptr = s1r_ptr + max_s1q
c        end do r = 0, s-1
         end do
            s1s_ptr = s1s_ptr + max_s1q
c        end do s = 1, max_s-1
         end do

c   o close the loop over the S2 ket irreps
c     end if (max_s.gt.1)
      end if
c     end do irp_s = 1, nirrep
      end do

c ----------------------------------------------------------------------

c     else if (irp_rs.ne.1) then
      else

c   o loop over the S2 ket irreps
      do irp_s = 1, nirrep
         max_s = s1_ket_pop(irp_s)
      if (max_s.ne.0) then
         irp_r = dirprd(irp_s,irp_rs)
         if (irp_r.lt.irp_s) then
            max_r = s1_ket_pop(irp_r)
            do_irp_s = (max_r.ne.0)
         else
            do_irp_s = .false.
         end if
      if (do_irp_s) then

c      o loop over the r<s orbital pairs (as offsets)
         to_r = max_r-1
         do s = 0, max_s-1
         do r = 0, to_r

c         o do work if rs is above the low bound
            if (rs.gt.rs_lo) then

c            o loop over the S2 bra irreps
               do irp_q = 1, nirrep
                  max_q = s1_bra_pop(irp_q)
               if (max_q.ne.0) then
                  irp_p = dirprd(irp_q,irp_rs)
                  if (irp_p.lt.irp_q) then
                     max_p = s1_bra_pop(irp_p)
                     do_irp_q = (max_p.ne.0)
                  else
                     do_irp_q = .false.
                  end if
               if (do_irp_q) then

c ----------------------------------------------------------------------

c   o S2(pq,rs) += alpha*S1(p,r)*S1(q,s)
      if ((irp_p.eq.dirprd(irp_r,irp_s1)).and.
     &    (irp_q.eq.dirprd(irp_s,irp_s1))     ) then

c      o point to the first S1(p,r) and S1(q,s) elements
         s1r_ptr = 1 + off_s1(irp_r) + max_p*r
         s1s_ptr = 1 + off_s1(irp_s) + max_q*s

c      o multiply and increment
#ifdef _USE_XGEMM
         call xgemm('T','N',max_p,max_q,1,
     &              alpha,s1(s1r_ptr),1,
     &                    s1(s1s_ptr),1,
     &              1.d0, s2(s2_ndx), max_p)
         s2_ndx = s2_ndx + max_p*max_q
#else
         to_p = max_p-1
         do q = 0, max_q-1
            dTmp_s = alpha*s1(s1s_ptr+q)
            do p = 0, to_p
               s2(s2_ndx+p) = s2(s2_ndx+p) + dTmp_s*s1(s1r_ptr+p)
            end do
            s2_ndx = s2_ndx + max_p
         end do
#endif /* _USE_XGEMM */

c     else if ((P != R x S1) || (Q != S x S1)) then
      else

c   o S2(pq,rs) -= alpha*S1(q,r)*S1(p,s)
      if ((irp_q.eq.dirprd(irp_r,irp_s1)).and.
     &    (irp_p.eq.dirprd(irp_s,irp_s1))     ) then

c      o point to the first S1(q,r) and S1(p,s) elements
         s1r_ptr = 1 + off_s1(irp_r) + max_q*r
         s1s_ptr = 1 + off_s1(irp_s) + max_p*s

c      o multiply and increment
#ifdef _USE_XGEMM
         call xgemm('T','N',max_p,max_q,1,
     &              -alpha,s1(s1s_ptr),1,
     &                     s1(s1r_ptr),1,
     &              1.d0,  s2(s2_ndx), max_p)
         s2_ndx = s2_ndx + max_p*max_q
#else
         to_p = max_p-1
         do q = 0, max_q-1
            dTmp_r = -alpha*s1(s1r_ptr+q)
            do p = 0, to_p
               s2(s2_ndx+p) = s2(s2_ndx+p) + dTmp_r*s1(s1s_ptr+p)
            end do
            s2_ndx = s2_ndx + max_p
         end do
#endif /* _USE_XGEMM */

c     else if (increment S2 manually) then
      else

c      o point to the next set of S2 p<q pairs
         s2_ndx = s2_ndx + max_p*max_q

c   o close the symmetry conditional
c     end if ((Q == R x S1) && (P == S x S1))
      end if

c   o close the symmetry conditional
c     end if ((P == R x S1) && (Q == S x S1))
      end if

c ----------------------------------------------------------------------

c            o close the loop over the S2 bra irreps
c              end if (do_irp_q)
               end if
c              end if (max_q.ne.0)
               end if
c              end do irp_q = 1, nirrep
               end do

c         o close the low bound conditional
c           end if (rs.gt.rs_lo)
            end if

c         o increment the rs super-index and quit if necessary
            rs = rs + 1
            if (rs.eq.rs_hi) return

c      o close the loop over the r<s orbital pairs
c        end do r = 0, to_r
         end do
c        end do s = 0, max_s-1
         end do

c   o close the loop over the S2 ket irreps
c     end if (do_irp_s)
      end if
c     end if (max_s.ne.0)
      end if
c     end do irp_s = 1, nirrep
      end do

c ----------------------------------------------------------------------

c     end if (irp_rs.eq.1)
      end if

      return
c     end subroutine antisym_s2pack_s1pq_sqr
      end

