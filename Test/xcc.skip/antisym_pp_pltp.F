
c This routine antisymmetrizes lists of distributions of one symmetry
c from <pq| for all p and q to <pq|| for p<q ( = <pq| - <qp| ).
c Orbitals p and q must span the same space.

c INPUT
c    double precision gs(*) : the source distributions

c    integer pq_size     : the length of one <pq| distribution
c    integer orb(nirrep) : the orbital occupation vector
c    integer nlists      : the number of distributions to process
c    integer irp_dist    : the symmetry of the distributions
c    integer offsets(nirrep*nirrep) : the matrix/array of offsets
c                                     to address blocks of unpacked irreps
c                                     (for any irp_dist)

c OUTPUT
c    double precision hh(*) : the antisymmetrized distributions

c ----------------------------------------------------------------------

      subroutine antisym_pp_pltp(gs,pq_size,hh,
     &                           orb,nlists,irp_dist,offsets)
      implicit none

c ARGUMENT LIST
      double precision gs(*), hh(*)
      integer pq_size, orb(*), nlists, irp_dist, offsets(*)

c INTERNAL VARIABLES
      integer add_gs, add_hh
      integer list, ioff
      integer p, irp_p, max_p
      integer q, irp_q, max_q
      logical do_irp_q

c COMMON BLOCKS
#include "syminf.com" /* for dirprd and nirrep */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      list = 0
c   o assert the number of lists is non-negative
      if (nlists.lt.0) then
         print *, '@ANTISYM_PP_PLTP: Assertion failed.'
         print *, '   nlists = ',nlists
         list = -1
      end if
c   o assert the distribution size is non-negative
      if (pq_size.lt.0) then
         print *, '@ANTISYM_PP_PLTP: Assertion failed.'
         print *, '   pq_size = ',pq_size
         list = -1
      end if
c   o assert irp_dist is in [1,nirrep]
      if ((irp_dist.lt.1).or.(nirrep.lt.irp_dist)) then
         print *, '@ANTISYM_PP_PLTP: Assertion failed.'
         print *, '   irp_dist = ',irp_dist
         list = -1
      end if
      if (list.ne.0) then
         call xcc_errex
      end if
#endif /* _ASSERT */

      if (nlists.gt.0) then

      add_gs = 1
      add_hh = 1

c   o do the symmetric point group separately
      if (irp_dist.eq.1) then

         do list = 1, nlists

            do irp_q = 1, nirrep
               max_q = orb(irp_q)
               do_irp_q = (max_q.gt.1)
            if (do_irp_q) then

c            o copy <pq|
               ioff = 0
               do q = 1, max_q-1
                  do p = 0, q-1
                     hh(add_hh+ioff) = gs(add_gs+(max_q*q)+p)
                     ioff = ioff + 1
                  end do
               end do

cYAU Replaced above hh(add_hh) with hh(add_hh+ioff) making this unnecessary.
cYAU c            o reset the destination pointer
cYAU                add_hh = add_hh - rshift(max_q*(max_q-1),1)

c            o subtract <qp|
               do q = 1, max_q-1
                  do p = 0, q-1
                     hh(add_hh) =   hh(add_hh)
     &                            - gs(add_gs+(max_q*p)+q)
                     add_hh = add_hh + 1
                  end do
               end do

c           end if (do_irp_q)
            end if
               add_gs = add_gs + ( max_q * max_q )
c           end do irp_q = 1, nirrep
            end do

c        end do list = 1, nlists
         end do

c     else if (irp_dist.ne.1) then
      else

         do list = 1, nlists

            do irp_q = 1, nirrep
               irp_p = dirprd(irp_q,irp_dist)
               if (irp_p.lt.irp_q) then
                  max_p = orb(irp_p)
                  max_q = orb(irp_q)
                  do_irp_q = ((max_p.ne.0).and.(max_q.ne.0))
               else
                  do_irp_q = .false.
               end if
            if (do_irp_q) then

c            o copy <pq|
               ioff = offsets(nirrep*(irp_q-1)+irp_p)
               do p = 0, (max_p*max_q)-1
                  hh(add_hh+p) = gs(add_gs+ioff+p)
               end do

c            o subtract <qp|
               ioff = offsets(nirrep*(irp_p-1)+irp_q)
               do q = 0, max_q-1
                  do p = 0, max_p-1
                     hh(add_hh) =   hh(add_hh)
     &                            - gs(add_gs+ioff+q+(max_q*p))
                     add_hh = add_hh + 1
                  end do
               end do

c           end if (do_irp_q)
            end if
c           end do irp_q = 1, nirrep
            end do

c         o point to the next <pq| list
            add_gs = add_gs + pq_size

c        end do list = 1, nlists
         end do

c     end if (irp_dist.eq.1)
      end if

c     end if (nlists.gt.0)
      end if

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

