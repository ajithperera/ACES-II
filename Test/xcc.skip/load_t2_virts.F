
c This routine loads T2(IJ,PQ) or T2(Ij,Pq) for all P orbitals under a
c range of Q orbitals and a single IJ/Ij pair. Currently, this routine
c is too smart for its own good. The caller tells it what it wants
c and where to load the amplitudes, but this routine relies heavily
c on the contents of the xcc common block to figure out how to
c process the data independently of what the caller may tell it.

c INPUT
c    [yeah, right]

c ----------------------------------------------------------------------

c#define _DEBUG_LOAD_T2_VIRTS

      subroutine load_t2_virts(
     &              dist, dest,
     &              irp_p, max_p,               ispin1,
     &              irp_q, max_q, from_q, to_q, ispin2, off_t2pq
     &           )
      implicit none

c ARGUMENT LIST
      double precision dist(*), dest(*)
      integer irp_p, max_p, ispin1
      integer irp_q, max_q, from_q, to_q, ispin2, off_t2pq

c INTERNAL VARIABLES
      integer p, q
      integer itmp, iadd, ndx

c COMMON BLOCKS
#include "xcc.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      itmp = 0
c   o assert (0 <= from_q <= to_q < max_q)
      if ((max_q.le.to_q).or.(to_q.lt.from_q).or.(from_q.lt.0)) then
         print *, '@LOAD_T2_VIRTS: Assertion failed.'
         print *, '   from_q = ',from_q
         print *, '   to_q   = ',to_q
         print *, '   max_q  = ',max_q
         itmp = 1
      end if
c   o assert (0 <= max_p)
      if (max_p.lt.0) then
         print *, '@LOAD_T2_VIRTS: Assertion failed.'
         print *, '   max_p = ',max_p
         itmp = 1
      end if
c   o assert ispin1 and ispin2 are 1 or 2
      if ((ispin1.lt.1).or.(ispin1.gt.2).or.
     &    (ispin2.lt.1).or.(ispin2.gt.2)    ) then
         print *, '@LOAD_T2_VIRTS: Assertion failed.'
         print *, '   ispin1 = ',ispin1
         print *, '   ispin2 = ',ispin2
         itmp = 1
      end if
c   o assert irp_p and irp_q are in [1,nirrep]
      if ((irp_p.lt.1).or.(nirrep.lt.irp_p).or.
     &    (irp_q.lt.1).or.(nirrep.lt.irp_q)    ) then
         print *, '@LOAD_T2_VIRTS: Assertion failed.'
         print *, '   irp_p = ',irp_p
         print *, '   irp_q = ',irp_q
         itmp = 1
      end if
c   o assert something about off_t2pq?
      if (itmp.ne.0) call xcc_errex
#endif /* _ASSERT */

      if (max_p.lt.1) return

c ----------------------------------------------------------------------

      if (ispin1.eq.ispin2) then
         if (irp_p.eq.irp_q) then
c            itmp = rshift(irp_q*(irp_q-1),1) + irp_q
            itmp = rshift(irp_q*(irp_q+1),1)
            iadd = 1 + off_pltp(itmp,off_t2pq)
            call unpack_pltp_cols(
     &              dist(iadd), dest,
     &              from_q, to_q, max_q-1
     &           )
c        else if (irp_p.ne.irp_q) then
         else
            if (irp_p.lt.irp_q) then
               itmp = rshift(irp_q*(irp_q-1),1) + irp_p
               iadd = 1 + off_pltp(itmp,off_t2pq) + ( from_q * max_p )
               itmp = max_p * ( 1 + to_q - from_q )
               do p = 0, itmp-1
                  dest(1+p) = dist(iadd+p)
               end do
            else
               itmp = rshift(irp_p*(irp_p-1),1) + irp_q
               itmp = 1 + off_pltp(itmp,off_t2pq) + from_q
               ndx = 1
               do q = from_q, to_q
                  iadd = itmp
                  do p = 0, max_p-1
                     dest(ndx+p) = -dist(iadd)
                     iadd = iadd + max_q
                  end do
                  ndx = ndx + max_p
                  itmp = itmp + 1
               end do
            end if
c        end if (irp_p.eq.irp_q)
         end if
c     else if (ispin1.ne.ispin2) then
      else
         if (ispin1.eq.1) then
            itmp = irp_p + nirrep*(irp_q-1)
            iadd = 1 + off_pq(itmp,h_Offsets_Ab) + ( from_q * max_p )
            itmp = max_p * ( 1 + to_q - from_q )
            do p = 0, itmp-1
               dest(1+p) = dist(iadd+p)
            end do
c        else if (ispin1.eq.2) then
         else
            itmp = irp_q + nirrep*(irp_p-1)
            itmp = 1 + off_pq(itmp,h_Offsets_Ab) + from_q
            ndx  = 1
            do q = from_q, to_q
               iadd = itmp
               do p = 0, max_p-1
                  dest(ndx+p) = -dist(iadd)
                  iadd = iadd + max_q
               end do
               ndx = ndx + max_p
               itmp = itmp + 1
            end do
c        end if (ispin1.eq.1)
         end if
c     end if (ispin1.eq.ispin2)
      end if

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

