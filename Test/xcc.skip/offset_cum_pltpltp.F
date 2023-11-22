
c This routine returns the cumulative offsets from element (1,1,1) needed
c to address any combination of symmetry-packed irreps of symmetry-packed
c orbitals of a 3-particle operator given a total symmetry irp_pqr.
c EXAMPLE: Starting from the first element in a list of T3 amplitudes
c          for particular i, j, and k orbitals, to reach the beginning
c          of the block:
c             irp_c = 3
c             irp_b = 3
c             irp_a = dirprd(irp_bc,irp_pqr)
c          you need to skip offsets(6) values provided pop()=vrt(1,ispin).
c The final value, offsets(1+nirrep*(nirrep+1)/2) contains the total number
c of P<Q<R elements in the list.

c INPUT
c    integer irp_pqr     : the total symmetry of the bra or ket
c    integer pop(nirrep) : the symmetry population vector

c OUTPUT
c    integer offsets(1+nirrep*(nirrep+1)/2) :
c                          the offsets for irp_p, irp_q, and irp_r
c                          in upper-triangular packed form
c                          addressable by irp_q+(irp_r*(irp_r-1)/2)

      subroutine offset_cum_pltpltp(irp_pqr,pop,offsets)
      implicit none

c ARGUMENT LIST
      integer irp_pqr, pop(*), offsets(*)

c INTERNAL VARIABLES
      integer irp_p, irp_q, irp_r, irp_qr
      integer max_p, max_q, max_r
      integer itmp, iadd

c COMMON BLOCKS
#include "syminf.com" /* sets nirrep and dirprd(,) */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      offsets(1) = 0
      iadd = 1
      do irp_r = 1, nirrep
         max_r = pop(irp_r)
         do irp_q = 1, irp_r
            irp_qr = dirprd(irp_q,irp_r)
            irp_p  = dirprd(irp_qr,irp_pqr)
            if (irp_p.eq.irp_q) then
               if (irp_q.eq.irp_r) then
                  itmp = ( max_r * (max_r-1) * (max_r-2) ) / 6
               else
                  max_q = pop(irp_q)
                  itmp  = ( max_r * max_q * (max_q-1) ) / 2
               end if
               offsets(iadd+1) = offsets(iadd) + itmp
            else
               if (irp_p.lt.irp_q) then
                  max_p = pop(irp_p)
                  if (irp_q.eq.irp_r) then
                     itmp = ( max_r * (max_r-1) * max_p ) / 2
                  else
                     max_q = pop(irp_q)
                     itmp  = ( max_r * max_q * max_p )
                  end if
                  offsets(iadd+1) = offsets(iadd) + itmp
               else
                  offsets(iadd+1) = offsets(iadd)
               end if
c           end if (irp_p.eq.irp_q)
            end if
            iadd = iadd + 1
c        end do irp_q = 1, irp_r
         end do
c     end do irp_r = 1, nirrep
      end do

c      print *, 'The P<Q<R offset vector is:'
c      print *, (offsets(itmp),itmp=1,1+(nirrep*(nirrep+1)/2))

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

