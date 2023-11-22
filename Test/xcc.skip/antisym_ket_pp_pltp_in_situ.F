
c This routine antisymmetrizes whole distributions of distributions
c for one symmetry from <*|pp> for all p (unpacked) to <*||pp> for p<p.
c The transformation is done in place destroying the source.

c EXAMPLE
c    For RHF: <Ka|Ij>;Ka,Ij -> <KA||IJ>;KA,I<J

c INPUT
c    double precision gs(*) : the source distributions

c    integer len_bra     : the length of one bra distribution
c    integer orb(nirrep) : the orbital occupation vector
c    integer irp_ket     : the symmetry of the ket distributions
c    integer offsets(nirrep*nirrep) : the matrix/array of offsets
c                                     to address blocks of unpacked irreps
c                                     (for any irp_ket)

c OUTPUT
c    double precision gs(*) : the antisymmetrized distributions

c ----------------------------------------------------------------------

      subroutine antisym_ket_pp_pltp_in_situ(gs,len_bra,
     &                                       orb,irp_ket,
     &                                       offsets)
      implicit none

c ARGUMENT LIST
      double precision gs(*)
      integer len_bra, orb(*), irp_ket, offsets(*)

c INTERNAL VARIABLES
      integer add_gs, add_hh
      integer p, irp_p, max_p, start_p, add_p, inc_p
      integer q, irp_q, max_q, ndx
      integer b
      logical do_irp_q

c COMMON BLOCKS
#include "syminf.com" /* for dirprd and nirrep */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      ndx = 0
c   o assert the distribution size is non-negative
      if (len_bra.lt.0) then
         print *, '@ANTISYM_KET_PP_PLTP_IN_SITU: Assertion failed.'
         print *, '   len_bra = ',len_bra
         ndx = -1
      end if
c   o assert irp_ket is in [1,nirrep]
      if ((irp_ket.lt.1).or.(nirrep.lt.irp_ket)) then
         print *, '@ANTISYM_KET_PP_PLTP_IN_SITU: Assertion failed.'
         print *, '   irp_ket = ',irp_ket
         ndx = -1
      end if
      if (ndx.ne.0) then
         call xcc_errex
      end if
#endif /* _ASSERT */

      if (len_bra.gt.0) then

      if (irp_ket.ne.1) then

c      o antisymmetrize all the p<q elements
         do irp_q = 1, nirrep
            irp_p = dirprd(irp_q,irp_ket)
            if (irp_p.lt.irp_q) then
               max_p = orb(irp_p)
               max_q = orb(irp_q)
               do_irp_q = ((max_p.ne.0).and.(max_q.ne.0))
            else
               do_irp_q = .false.
            end if
         if (do_irp_q) then

c         o |pq> = |pq> - |qp>
            start_p = 1 + ( len_bra * offsets(nirrep*(irp_p-1)+irp_q) )
            ndx     = 1 + ( len_bra * offsets(nirrep*(irp_q-1)+irp_p) )
            inc_p   = len_bra * max_q
            do q = 0, max_q-1
               add_p = start_p
               do p = 0, max_p-1
                  do b = 0, len_bra-1
                     gs(ndx+b) = gs(ndx+b) - gs(add_p+b)
                  end do
                  ndx = ndx + len_bra
                  add_p = add_p + inc_p
               end do
               start_p = start_p + len_bra
            end do

c        end if (do_irp_q)
         end if
c        end do irp_q = 1, nirrep
         end do

c      o pack the distributions
         add_hh = 1
         do irp_q = 1, nirrep
            irp_p = dirprd(irp_q,irp_ket)
            if (irp_p.lt.irp_q) then
               max_p = orb(irp_p)
               max_q = orb(irp_q)
               do_irp_q = ((max_p.ne.0).and.(max_q.ne.0))
            else
               do_irp_q = .false.
            end if
         if (do_irp_q) then

            add_gs = 1 + ( len_bra * offsets(nirrep*(irp_q-1)+irp_p) )
            ndx = len_bra * max_p * max_q
            do b = 0, ndx-1
               gs(add_hh+b) = gs(add_gs+b)
            end do
            add_hh = add_hh + ndx

c        end if (do_irp_q)
         end if
c        end do irp_q = 1, nirrep
         end do

c   o do the symmetric point group separately
c     else if (irp_ket.eq.1) then
      else

         add_gs = 1
         add_hh = 1
         do irp_q = 1, nirrep
            max_q = orb(irp_q)
            do_irp_q = (max_q.gt.1)
         if (do_irp_q) then

c         o |pq> = |pq> - |qp> (for p<q)
            start_p = add_gs + len_bra
            inc_p   = len_bra * max_q
            do q = 1, max_q-1
               ndx = add_gs + ( len_bra * max_q * q )
               add_p = start_p
               do p = 0, q-1
                  do b = 0, len_bra-1
                     gs(ndx+b) = gs(ndx+b) - gs(add_p+b)
                  end do
                  ndx = ndx + len_bra
                  add_p = add_p + inc_p
               end do
               start_p = start_p + len_bra
            end do

c         o pack the list
            inc_p = len_bra
            do q = 1, max_q-1
               ndx = add_gs + ( len_bra * max_q * q )
               do b = 0, inc_p-1
                  gs(add_hh+b) = gs(ndx+b)
               end do
               add_hh = add_hh + inc_p
               inc_p  = inc_p  + len_bra
            end do

c        end if (do_irp_q)
         end if
            add_gs = add_gs + ( len_bra * max_q * max_q )
c        end do irp_q = 1, nirrep
         end do

c     end if (irp_ket.ne.1)
      end if

c     end if (len_bra.gt.0)
      end if

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

