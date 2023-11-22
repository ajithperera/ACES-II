
c This routine loads T2(MI,PQ) or T2(Mi,Pq) for all P<Q or Pq
c under a given range of M and a single I/i. Currently, this routine
c is too smart for its own good. The caller tells it what it wants
c and where to load the amplitudes, but this routine relies heavily
c on the contents of the xcc common block to figure out how to
c process the data independently of what the caller may tell it.

c INPUT
c    [yeah, right]

c ----------------------------------------------------------------------

c#define _DEBUG_LOAD_T2_DISTS

      subroutine load_t2_dists(
     &              dcore, dest,
     &              irp_pq, max_pq,
     &              irp_m,  max_m,  from_m, to_m, ispin1,
     &              irp_i,  max_i,  i,            ispin2, off_t2mi
     &           )
      implicit none

c ARGUMENT LIST
      double precision dcore(*), dest(*)
      integer irp_pq, max_pq
      integer irp_m,  max_m, from_m, to_m, ispin1
      integer irp_i,  max_i, i, ispin2, off_t2mi

c INTERNAL VARIABLES
      integer pq, m, num_m
      integer ioff, itmp, iadd, ndx

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
c   o assert (0 <= from_m <= to_m < max_m)
      if ((max_m.le.to_m).or.(to_m.lt.from_m).or.(from_m.lt.0)) then
         print *, '@LOAD_T2_DISTS: Assertion failed.'
         print *, '   from_m = ',from_m
         print *, '   to_m   = ',to_m
         print *, '   max_m  = ',max_m
         itmp = 1
      end if
c   o assert (0 <= max_pq)
      if (max_pq.lt.0) then
         print *, '@LOAD_T2_DISTS: Assertion failed.'
         print *, '   max_pq = ',max_pq
         itmp = 1
      end if
c   o assert ispin1 and ispin2 are 1 or 2
      if ((ispin1.lt.1).or.(ispin1.gt.2).or.
     &    (ispin2.lt.1).or.(ispin2.gt.2)    ) then
         print *, '@LOAD_T2_DISTS: Assertion failed.'
         print *, '   ispin1 = ',ispin1
         print *, '   ispin2 = ',ispin2
         itmp = 1
      end if
c   o assert i is well-behaved
      if ((max_i.le.i).or.(i.lt.0)) then
         print *, '@LOAD_T2_DISTS: Assertion failed.'
         print *, '   i = ',i,', max_i = ',max_i
         itmp = 1
      end if
c   o assert irp_m, irp_i, and irp_pq are in [1,nirrep]
      if ((irp_m .lt.1).or.(nirrep.lt.irp_m ).or.
     &    (irp_i .lt.1).or.(nirrep.lt.irp_i ).or.
     &    (irp_pq.lt.1).or.(nirrep.lt.irp_pq)    ) then
         print *, '@LOAD_T2_DISTS: Assertion failed.'
         print *, '   irp_m  = ',irp_m
         print *, '   irp_i  = ',irp_i
         print *, '   irp_pq = ',irp_pq
         itmp = 1
      end if
c   o assert irp_pq is irp_m x irp_i
      if (irp_pq.ne.dirprd(irp_m,irp_i)) then
         print *, '@LOAD_T2_DISTS: Assertion failed.'
         print *, '   irp_pq = ',irp_pq
         print *, '   irp_m  = ',irp_m
         print *, '   irp_i  = ',irp_i
         print *, '   irp_m x irp_i = ',dirprd(irp_m,irp_i)
         itmp = 1
      end if
c   o assert something about off_t2mi?
      if (itmp.ne.0) call xcc_errex
#endif /* _ASSERT */

      if (max_pq.lt.1) return

c ----------------------------------------------------------------------

c   o get T2(MI,PQ) distributions (all P<Q for all M)
      if (ispin1.eq.ispin2) then
         if (irp_m.eq.irp_i) then
            ndx  = 1
c            itmp = rshift(irp_i*(irp_i-1),1) + irp_i
            itmp = rshift(irp_i*(irp_i+1),1)
            ioff = off_pltp(itmp,off_t2mi) + rshift(i*(i-1),1)
            iadd = ndx_t2(irp_pq,ispin1) + ( ioff * max_pq )
            if (i.ne.0) then
               if (from_m.lt.i) then
c               o load all M<I
                  iadd = iadd + ( from_m * max_pq )
                  itmp = ( min(i,1+to_m) - from_m ) * max_pq
                  do pq = 0, itmp-1
                     dest(ndx+pq) = dcore(iadd+pq)
                  end do
                  ndx  = ndx  + itmp
                  iadd = iadd + itmp + ( i * max_pq )
               else
                  iadd = iadd + ( lshift(i,1) * max_pq )
               end if
            end if
            if ((from_m.lt.(1+i)).and.(i.lt.(1+to_m))) then
c            o T2(II,PQ) is 0
               do pq = 0, max_pq-1
                  dest(ndx+pq) = 0.0d0
               end do
               ndx = ndx + max_pq
            end if
            if ((i.lt.to_m).and.((1+i).ne.max_i)) then
c            o load all M>I
               itmp = from_m-(1+i)
               do m = 1, min(max_i-1,to_m)-i
                  if (m.gt.itmp) then
                     do pq = 0, max_pq-1
                        dest(ndx+pq) = -dcore(iadd+pq)
                     end do
                     ndx = ndx + max_pq
                  end if
                  iadd = iadd + ( m * max_pq ) + ( i * max_pq )
               end do
            end if
c        else if (irp_m.ne.irp_i) then
         else
            num_m = 1 + to_m - from_m
            if (irp_m.lt.irp_i) then
               itmp = rshift(irp_i*(irp_i-1),1) + irp_m
               ioff = off_pltp(itmp,off_t2mi) + ( i * max_m ) + from_m
               iadd = ndx_t2(irp_pq,ispin1) + ( ioff * max_pq )
               itmp = num_m * max_pq
               do pq = 0, itmp-1
                  dest(1+pq) = dcore(iadd+pq)
               end do
            else
               itmp = rshift(irp_m*(irp_m-1),1) + irp_i
               ioff = off_pltp(itmp,off_t2mi) + ( from_m * max_i ) + i
               iadd = ndx_t2(irp_pq,ispin1) + ( ioff * max_pq )
               itmp = max_i * max_pq
               ndx  = 1
               do m = 1, num_m
                  do pq = 0, max_pq-1
                     dest(ndx+pq) = -dcore(iadd+pq)
                  end do
                  ndx  = ndx  + max_pq
                  iadd = iadd + itmp
               end do
            end if
c        end if (irp_m.eq.irp_i)
         end if
c     else if (ispin1.ne.ispin2) then
      else
         num_m = 1 + to_m - from_m
         if (ispin1.eq.1) then
            itmp = irp_m + nirrep*(irp_i-1)
            ioff = off_pq(itmp,h_Offsets_Ij) + ( i * max_m ) + from_m
            iadd = ndx_t2(irp_pq,3) + ( ioff * max_pq )
            itmp = num_m * max_pq
            do pq = 0, itmp-1
               dest(1+pq) = dcore(iadd+pq)
            end do
c        else if (ispin1.eq.2) then
         else
            itmp = irp_i + nirrep*(irp_m-1)
            ioff = off_pq(itmp,h_Offsets_Ij) + ( from_m * max_i ) + i
            iadd = ndx_t2(irp_pq,3) + ( ioff * max_pq )
            itmp = max_i * max_pq
            ndx  = 1
            do m = 1, num_m
               do pq = 0, max_pq-1
                  dest(ndx+pq) = dcore(iadd+pq)
               end do
               ndx  = ndx  + max_pq
               iadd = iadd + itmp
            end do
            call transp_dist(dest,vrt(1,1),vrt(1,2),irp_pq,num_m,
     &                       .false.)
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

