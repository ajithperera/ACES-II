
c This routine returns the cumulative offsets from element (1,1,1) needed
c to address any combination of symmetry-packed irreps of symmetry-packed
c orbitals of a 3-particle operator given a total symmetry irp_pqr.
c EXAMPLE: Starting from the first element in a list of T3 amplitudes
c          for particular I, J, and k orbitals, to reach the beginning
c          of the block:
c             irp_c = 3
c             irp_B = 3
c             irp_A = dirprd(irp_bc,irp_pqr)
c          you need to skip offsets(nirrep*2+3) values provided
c          pop1()=vrt(1,1) and pop2()=vrt(1,2).
c The final value, offsets(1+nirrep*nirrep) contains the total number of
c P<Q,r elements in the list.

c INPUT
c    integer irp_pqr      : the total symmetry of the bra or ket
c    integer pop1(nirrep) : the symmetry population vector of spin 1
c    integer pop2(nirrep) : the symmetry population vector of spin 2

c OUTPUT
c    integer offsets(1+nirrep*nirrep) :
c                           the offsets for irp_P, irp_Q, and irp_r
c                           addressable by nirrep*(irp_r-1)+irp_q

      subroutine offset_cum_pltpq(irp_pqr,pop1,pop2,offsets)
      implicit none

c ARGUMENT LIST
      integer irp_pqr, pop1(*), pop2(*), offsets(*)

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
         max_r = pop2(irp_r)
         do irp_q = 1, nirrep
            irp_qr = dirprd(irp_q,irp_r)
            irp_p  = dirprd(irp_qr,irp_pqr)
            if (irp_p.eq.irp_q) then
               max_q = pop1(irp_q)
               itmp  = ( max_r * max_q * (max_q-1) ) / 2
               offsets(iadd+1) = offsets(iadd) + itmp
            else
               if (irp_p.lt.irp_q) then
                  max_p = pop1(irp_p)
                  max_q = pop1(irp_q)
                  itmp  = ( max_r * max_q * max_p )
                  offsets(iadd+1) = offsets(iadd) + itmp
               else
                  offsets(iadd+1) = offsets(iadd)
               end if
c           end if (irp_p.eq.irp_q)
            end if
            iadd = iadd + 1
c        end do irp_q = 1, nirrep
         end do
c     end do irp_r = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

