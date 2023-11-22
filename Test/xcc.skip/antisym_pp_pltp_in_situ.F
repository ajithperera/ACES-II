
c This routine antisymmetrizes lists of distributions of one symmetry
c from <pq| for all p and q to <pq|| for p<q ( = <pq| - <qp| ).
c The transformation is done in place destroying the source.

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
c    double precision gs(*) : the antisymmetrized distributions

c ----------------------------------------------------------------------

      subroutine antisym_pp_pltp_in_situ(gs,pq_size,
     &                                   orb,nlists,irp_dist,offsets)
      implicit none

c ARGUMENT LIST
      double precision gs(*)
      integer pq_size, orb(*), nlists, irp_dist, offsets(*)

c INTERNAL VARIABLES
      integer add_gs, add_hh
      integer list, ioff_p, ioff_q
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
         print *, '@ANTISYM_PP_PLTP_IN_SITU: Assertion failed.'
         print *, '   nlists = ',nlists
         list = -1
      end if
c   o assert the distribution size is non-negative
      if (pq_size.lt.0) then
         print *, '@ANTISYM_PP_PLTP_IN_SITU: Assertion failed.'
         print *, '   pq_size = ',pq_size
         list = -1
      end if
c   o assert irp_dist is in [1,nirrep]
      if ((irp_dist.lt.1).or.(nirrep.lt.irp_dist)) then
         print *, '@ANTISYM_PP_PLTP_IN_SITU: Assertion failed.'
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

c            o <pq| = <pq| - <qp| (for p<q)
               do q = 1, max_q-1
                  do p = 0, q-1
                     gs(add_gs+(max_q*q)+p) =
     &               gs(add_gs+(max_q*q)+p)
     &             - gs(add_gs+q+(max_q*p))
                  end do
               end do

c            o pack the list into the proper address
               do q = 1, max_q-1
                  do p = 0, q-1
                     gs(add_hh) = gs(add_gs+(max_q*q)+p)
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

c         o antisymmetrize all the p<q elements for this list
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

c            o <pq| = <pq| - <qp|
               ioff_q = offsets(nirrep*(irp_q-1)+irp_p)
               ioff_p = offsets(nirrep*(irp_p-1)+irp_q)
               do q = 0, max_q-1
                  do p = 0, max_p-1
                     gs(add_gs+ioff_q+(max_p*q)+p) =
     &               gs(add_gs+ioff_q+(max_p*q)+p)
     &             - gs(add_gs+ioff_p+q+(max_q*p))
                  end do
               end do

c           end if (do_irp_q)
            end if
c           end do irp_q = 1, nirrep
            end do

c         o pack the list into the proper address
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

               ioff_q = offsets(nirrep*(irp_q-1)+irp_p)
               ioff_p = max_p * max_q
               do p = 0, ioff_p-1
                  gs(add_hh+p) = gs(add_gs+ioff_q+p)
               end do
               add_hh = add_hh + ioff_p

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

