
c This routine creates the product of two 1-particle operators. The 2-particle
c result is scaled and added into the 2-particle destination. The symmetry-
c forbidden elements of the product are not touched.

c The core of the subroutine performs the following:
c    S(pq,rs) += alpha * [ S(p,r)*S(q,s) ]

c INPUT
c int irp_pq   : the total bra symmetry of the S2 operator
c int irp_rs   : the total ket symmetry of the S2 operator
c int rs_start : the rs super-index to start on (think getlst columns)
c int num_rs   : the number of rs pairs to update
c double alpha : the S1(P,R)*S1(Q,S) scaling factor
c double s1pr(*) : the S1(P,R) matrix elements for all ket irreps
c double s1qs(*) : the S1(Q,S) matrix elements for all ket irreps
c int irp_s1pr : the total symmetry of the S1(P,R) operator
c int irp_s1qs : the total symmetry of the S1(Q,S) operator
c int pop_p(*) : the population vector for the S1(P,R) bra elements
c int pop_r(*) : the population vector for the S1(P,R) ket elements
c int pop_q(*) : the population vector for the S1(Q,S) bra elements
c int pop_s(*) : the population vector for the S1(Q,S) ket elements
c int off_s1pr(*) : the vector of offsets to reach each S1(P,R) ket irrep
c int off_s1qs(*) : the vector of offsets to reach each S1(P,R) ket irrep

c INPUT/OUTPUT
c double s2(*) : the array to update with S1(P,R)*S1(Q,S)

c#define _DEBUG_S2PQRS_PE_S1PR_S1QS
c#define _USE_XGEMM

      subroutine s2pqrs_pe_s1pr_s1qs(
     &                    s2,irp_pq,irp_rs,rs_start,num_rs,
     &              alpha,s1pr,irp_s1pr,pop_p,pop_r,off_s1pr,
     &                    s1qs,irp_s1qs,pop_q,pop_s,off_s1qs
     &           )
      implicit none

c ARGUMENT LIST
      double precision alpha, s1pr(*), s1qs(*), s2(*)
      integer irp_pq, irp_rs, rs_start, num_rs
      integer irp_s1pr, pop_p(*), pop_r(*), off_s1pr(*)
      integer irp_s1qs, pop_q(*), pop_s(*), off_s1qs(*)

c INTERNAL VARIABLES
      double precision dTmp_s
      integer s1pr_ptr2, s1pr_ptr, s1qs_ptr, s2_ndx
      integer rs_lo, rs, rs_hi
      integer irp_s1p, max_s1p
      integer irp_s1q, max_s1q
      integer irp_r, r, max_r, to_r
      integer irp_s, s, max_s
      integer irp_p, p, max_p, to_p
      integer irp_q, q, max_q
      integer iTmp

c COMMON BLOCKS
#include "syminf.com" /* for nirrep and dirprd(,) */

c ----------------------------------------------------------------------

#ifdef _ASSERT
      iTmp = 0
c   o assert irp_pq, irp_rs, irp_s1pr, and irp_s1qs are in [1,nirrep]
      if ((irp_pq  .lt.1).or.(nirrep.lt.irp_pq  ).or.
     &    (irp_rs  .lt.1).or.(nirrep.lt.irp_rs  ).or.
     &    (irp_s1pr.lt.1).or.(nirrep.lt.irp_s1pr).or.
     &    (irp_s1qs.lt.1).or.(nirrep.lt.irp_s1qs)    ) then
         print *, '@S2PQRS_PE_S1PR_S1QS: Assertion failed.'
         print *, '   irp_pq   = ',irp_pq
         print *, '   irp_rs   = ',irp_rs
         print *, '   irp_s1pr = ',irp_s1pr
         print *, '   irp_s1qs = ',irp_s1qs
         iTmp = 1
      end if
c   o assert total symmetries are allowed
      if (dirprd(irp_pq,irp_rs).ne.dirprd(irp_s1pr,irp_s1qs)) then
         print *, '@S2PQRS_PE_S1PR_S1QS: Assertion failed.'
         print *, '   irp_pq   = ',irp_pq
         print *, '   irp_rs   = ',irp_rs
         print *, '   irp_s1pr = ',irp_s1pr
         print *, '   irp_s1qs = ',irp_s1qs
         iTmp = 1
      end if
c   o assert rs_start and num_rs are properly bound
      if ((rs_start.lt.1).or.(num_rs.lt.0)) then
         print *, '@S2PQRS_PE_S1PR_S1QS: Assertion failed.'
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

c   o loop over the S2 ket irreps
      do irp_s = 1, nirrep
         max_s = pop_s(irp_s)
      if (max_s.ne.0) then
         irp_r = dirprd(irp_s,irp_rs)
         max_r = pop_r(irp_r)
      if (max_r.ne.0) then

c      o for these ket irreps
c         - cache the direct product look-up for the symmetry-allowed bra irreps
c         - cache the populations for the bra irreps
c         - start the S1 column pointers
         irp_s1p = dirprd(irp_r,irp_s1pr)
         irp_s1q = dirprd(irp_s,irp_s1qs)
         max_s1p = pop_p(irp_s1p)
         max_s1q = pop_q(irp_s1q)
         s1pr_ptr2 = 1 + off_s1pr(irp_r)
         s1qs_ptr  = 1 + off_s1qs(irp_s)

c      o loop over the rs orbital pairs (as offsets)
         to_r = max_r-1
         do s = 0, max_s-1
c         o reset the S1(P,R) pointer
            s1pr_ptr = s1pr_ptr2
         do r = 0, to_r

c         o do work if rs is above the low bound
            if (rs.gt.rs_lo) then

c            o loop over the S2 bra irreps
               do irp_q = 1, nirrep
                  max_q = pop_q(irp_q)
               if (max_q.ne.0) then
                  irp_p = dirprd(irp_q,irp_pq)
                  max_p = pop_p(irp_p)
               if (max_p.ne.0) then

c ----------------------------------------------------------------------

c   o process symmetry-allowed elements
      if ((irp_p.eq.irp_s1p).and.
     &    (irp_q.eq.irp_s1q)     ) then

c      o multiply and increment
#ifdef _USE_XGEMM
         call xgemm('T','N',max_p,max_q,1,
     &              alpha,s1(s1pr_ptr),1,
     &                    s1(s1qs_ptr),1,
     &              1.d0, s2(s2_ndx),  max_p)
         s2_ndx = s2_ndx + max_p*max_q
#else
         to_p = max_p-1
         do q = 0, max_q-1
            dTmp_s = alpha*s1qs(s1qs_ptr+q)
            do p = 0, to_p
               s2(s2_ndx+p) = s2(s2_ndx+p) + dTmp_s*s1pr(s1pr_ptr+p)
            end do
            s2_ndx = s2_ndx + max_p
         end do
#endif /* _USE_XGEMM */

c     else if ((P != R x S1PR) && (Q != S x S1QS)) then
      else

c      o point to the next set of S2 pq pairs
         s2_ndx = s2_ndx + max_p*max_q

c   o close the symmetry conditional
c     end if ((P == R x S1PR) && (Q == S x S1QS))
      end if

c ----------------------------------------------------------------------

c            o close the loop over the S2 bra irreps
c              end if (max_p.ne.0)
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

c      o increment the S1 pointers and close the loop over the rs orbital pairs
            s1pr_ptr = s1pr_ptr + max_s1p
c        end do r = 0, to_r
         end do
            s1qs_ptr = s1qs_ptr + max_s1q
c        end do s = 0, max_s-1
         end do

c   o close the loop over the S2 ket irreps
c     end if (max_r.ne.0)
      end if
c     end if (max_s.ne.0)
      end if
c     end do irp_s = 1, nirrep
      end do

c ----------------------------------------------------------------------

      return
c     end subroutine s2pqrs_pe_s1pr_s1qs
      end

