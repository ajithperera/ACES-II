
c This routine returns the size of all P<Q<R distributions
c given a single population vector.

c INPUT
c    integer orb(nirrep) : a population vector of irreps

c OUTPUT
c    integer dpd(nirrep) : the P<Q<R distributions per irrep

      subroutine irpdpd_pltpltp(orb,dpd)
      implicit none

c ARGUMENT LIST
      integer orb(*), dpd(*)

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
            max_r = orb(irp_r)
            do irp_q = 1, irp_r
               irp_qr = dirprd(irp_q,irp_r)
               irp_p  = dirprd(irp_qr,irp_pqr)
               if (irp_p.eq.irp_q) then
                  if (irp_q.eq.irp_r) then
                     itmp = ( max_r * (max_r-1) * (max_r-2) ) / 6
                  else
                     max_q = orb(irp_q)
                     itmp  = ( max_r * max_q * (max_q-1) ) / 2
                  end if
                  dpd(irp_pqr) = dpd(irp_pqr) + itmp
               else
                  if (irp_p.lt.irp_q) then
                     max_p = orb(irp_p)
                     if (irp_q.eq.irp_r) then
                        itmp = ( max_r * (max_r-1) * max_p ) / 2
                     else
                        max_q = orb(irp_q)
                        itmp  = ( max_r * max_q * max_p )
                     end if
                     dpd(irp_pqr) = dpd(irp_pqr) + itmp
                  end if
c              end if (irp_p.eq.irp_q)
               end if
c           end do irp_q = 1, irp_r
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

