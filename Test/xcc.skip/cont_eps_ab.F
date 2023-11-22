
c This routine divides one operator distribution by epsilon. Specifically,
c the distribution must be stored all AB (or Ab). Additionally, the
c quantity O(AB,?)**2/eps_?_AB is returned. For example, if O is the
c 2-electron antisymmetrized integral, then the energy returned is a
c contribution to the MBPT(2) energy.

c INPUT:
c    double precision dist(*) : the distribution to modify

c    integer irp_ab : the symmetry of the distribution
c    integer vrt_a(*) : the population vector of the A virtual irreps
c    double precision eval_a(*)     : the virtual eigenvalues of A
c    integer          add_eval_a(*) : the starting addresses by irrep
c                                     of virtual eigenvalues in eval_a(*)
c    integer vrt_b(*) : the population vector of the B virtual irreps
c    double precision eval_b(*)     : the virtual eigenvalues of B
c    integer          add_eval_b(*) : the starting addresses by irrep
c                                     of virtual eigenvalues in eval_b(*)

c    double precision eps_ij : the positive part of the epsilon denominator

c OUTPUT:
c    double precision energy : the O**2/eps contribution (NOT ZEROED)

      subroutine cont_eps_ab(
     &              dist,
     &              irp_ab,vrt_a,eval_a,add_eval_a,
     &                     vrt_b,eval_b,add_eval_b,
     &              eps_ij,
     &              energy
     &           )
      implicit none

c ARGUMENT LIST
      double precision dist(*)
      integer irp_ab, vrt_a(*), vrt_b(*)
      double precision eval_a(*), eval_b(*)
      integer add_eval_a(*), add_eval_b(*)
      double precision eps_ij
      double precision energy

c INTERNAL VARIABLES
      integer a, irp_a, max_a
      integer b, irp_b, max_b
      logical do_irp_b
      integer ndx
      double precision eps_ij_b
      double precision eps_ij_ab
      double precision tmp_energy

c COMMON BLOCKS
#include "syminf.com" /* for nirrep and dirprd */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      ndx = 0
c   o assert irp_ab is in [1,nirrep]
      if ((irp_ab.lt.1).or.(nirrep.lt.irp_ab)) then
         print *, '@CONT_EPS_AB: Assertion failed.'
         print *, '   irp_ab = ',irp_ab
         ndx = -1
      end if
      if (ndx.ne.0) then
         call xcc_errex
      end if
#endif /* _ASSERT */

      tmp_energy = 0.0d0
      ndx        = 1

c   o loop over irreps of virtuals
      do irp_b = 1, nirrep
         irp_a = dirprd(irp_b,irp_ab)
         max_a = vrt_a(irp_a)
         max_b = vrt_b(irp_b)
         do_irp_b = ((max_b.ne.0).and.(max_a.ne.0))
      if (do_irp_b) then

c      o loop over orbitals AB
         do b = 0, max_b-1
            eps_ij_b  = eps_ij   - eval_b(add_eval_b(irp_b) + b)
         do a = 0, max_a-1
            eps_ij_ab = eps_ij_b - eval_a(add_eval_a(irp_a) + a)

c ----------------------------------------------------------------------

      dist(ndx) = dist(ndx) / eps_ij_ab
      tmp_energy = tmp_energy + ( dist(ndx) * dist(ndx) * eps_ij_ab )

c ----------------------------------------------------------------------

c      o increment and close the loop over orbitals AB
            ndx = ndx + 1
c        end do a = 0, max_a-1
         end do
c        end do b = 0, max_b-1
         end do

c   o close the loop over irreps of virtuals
c     end if (do_irp_b)
      end if
c     end do irp_b = 1, nirrep
      end do

      energy = energy + tmp_energy

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

