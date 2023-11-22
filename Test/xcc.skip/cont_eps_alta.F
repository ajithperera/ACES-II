
c This routine divides one operator distribution by epsilon. Specifically,
c the distribution must be stored all A<B where A and B span the same virtual
c space. Additionally, the quantity O(A<B,?)**2/eps_?_ab is returned.
c For example, if O is the 2-electron antisymmetrized integral, then the
c energy returned is a contribution to the MBPT(2) energy.

c INPUT:
c    double precision dist(*) : the distribution to modify

c    integer irp_ab : the symmetry of the distribution
c    integer vrt_a(*) : the population vector of the virtual irreps
c    double precision eval_a(*)     : the virtual eigenvalues
c    integer          add_eval_a(*) : the starting addresses by irrep
c                                     of virtual eigenvalues in eval_a(*)

c    double precision eps_ij  : the positive part of the epsilon denominator

c OUTPUT:
c    double precision energy : the O**2/eps contribution (NOT ZEROED)

c ----------------------------------------------------------------------

      subroutine cont_eps_alta(
     &              dist,
     &              irp_ab,vrt_a,eval_a,add_eval_a,
     &              eps_ij,
     &              energy
     &           )
      implicit none

c ARGUMENT LIST
      double precision dist(*)
      integer irp_ab, vrt_a(*)
      double precision eval_a(*)
      integer add_eval_a(*)
      double precision eps_ij
      double precision energy

c INTERNAL VARIABLES
      integer a, irp_a, max_a,         to_a
      integer b, irp_b, max_b, from_b
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
         print *, '@CONT_EPS_ALTA: Assertion failed.'
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
         if (irp_a.lt.irp_b) then
            max_a = vrt_a(irp_a)
            max_b = vrt_a(irp_b)
            do_irp_b = ((max_b.ne.0).and.(max_a.ne.0))
         else
            if (irp_a.gt.irp_b) then
               do_irp_b = .false.
            else
               max_a = vrt_a(irp_a)
               max_b = max_a
               do_irp_b = (max_b.gt.1)
            end if
         end if
      if (do_irp_b) then

c      o set orbital boundary conditions (as offsets)
         to_a = max_a-1
         if (irp_a.lt.irp_b) then
            from_b = 0
         else
            from_b = 1
         end if

c      o loop over orbitals A<B
         do b = from_b, max_b-1
            eps_ij_b  = eps_ij   - eval_a(add_eval_a(irp_b) + b)
            if (irp_a.eq.irp_b) to_a = b-1
         do a = 0, to_a
            eps_ij_ab = eps_ij_b - eval_a(add_eval_a(irp_a) + a)

c ----------------------------------------------------------------------

      dist(ndx) = dist(ndx) / eps_ij_ab
      tmp_energy = tmp_energy + ( dist(ndx) * dist(ndx) * eps_ij_ab )

c ----------------------------------------------------------------------

c      o increment and close the loop over orbitals A<B
            ndx = ndx + 1
c        end do a = 0, to_a
         end do
c        end do b = from_b, max_b-1
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

