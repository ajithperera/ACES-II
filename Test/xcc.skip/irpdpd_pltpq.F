
c This routine returns the size of all P<Q,r distributions
c given two population vectors.

c INPUT
c    integer orb1(nirrep) : a primary   population vector of irreps
c    integer orb2(nirrep) : a secondary population vector of irreps

c OUTPUT
c    integer dpd(nirrep) : the P<Q,r distributions per irrep

      subroutine irpdpd_pltpq(orb1,orb2,dpd)
      implicit none

c ARGUMENT LIST
      integer orb1(*), orb2(*), dpd(*)

c INTERNAL VARIABLES
      integer irp_p, irp_q, irp_r, irp_qr, irp_pqr
      integer max_p, max_q, max_r
      integer itmp

c COMMON BLOCKS
#include "syminf.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

      do irp_pqr = 1, nirrep
         dpd(irp_pqr) = 0
         do irp_r = 1, nirrep
            max_r = orb2(irp_r)
            do irp_q = 1, nirrep
               irp_qr = dirprd(irp_q,irp_r)
               irp_p  = dirprd(irp_qr,irp_pqr)
               if (irp_p.eq.irp_q) then
                  max_q = orb1(irp_q)
                  itmp  = ( max_r * max_q * (max_q-1) ) / 2
                  dpd(irp_pqr) = dpd(irp_pqr) + itmp
               else
                  if (irp_p.lt.irp_q) then
                     max_p = orb1(irp_p)
                     max_q = orb1(irp_q)
                     itmp  = ( max_r * max_q * max_p )
                     dpd(irp_pqr) = dpd(irp_pqr) + itmp
                  end if
c              end if (irp_p.eq.irp_q)
               end if
c           end do irp_q = 1, nirrep
            end do
c        end do irp_r = 1, nirrep
         end do
c     end do irp_pqr = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

