
c This routine divides one operator symmetry block by epsilon. Specifically,
c the distribution must be stored all AB (or Ab) for all IJ (or Ij).
c Additionally, the quantity O**2/eps is returned. For example, if O is the
c 2-electron antisymmetrized integral, then the energy returned is a
c contribution to the MBPT(2) energy.

c INPUT:
c    double precision dist(*) : the distributions to modify

c    integer irp_ab : the symmetry of the virtual distributions
c    integer vrt_a(*) : the population vector of the A virtual irreps
c    double precision eval_a(*)     : the virtual eigenvalues of A
c    integer          add_eval_a(*) : the starting addresses by irrep
c                                     of virtual eigenvalues in eval_a(*)
c    integer vrt_b(*) : the population vector of the B virtual irreps
c    double precision eval_b(*)     : the virtual eigenvalues of B
c    integer          add_eval_b(*) : the starting addresses by irrep
c                                     of virtual eigenvalues in eval_b(*)
c    integer dim_dist : the length of one AB distribution

c    integer irp_ij : the symmetry of the occupied distributions
c    integer occ_i(*) : the population vector of the I occupied irreps
c    double precision eval_i(*)     : the occupied eigenvalues of I
c    integer          add_eval_i(*) : the starting addresses by irrep
c                                     of occupied eigenvalues in eval_i(*)
c    integer occ_j(*) : the population vector of the J occupied irreps
c    double precision eval_j(*)     : the occupied eigenvalues of J
c    integer          add_eval_j(*) : the starting addresses by irrep
c                                     of occupied eigenvalues in eval_j(*)

c OUTPUT:
c    double precision energy : the O**2/eps contribution (NOT ZEROED)

      subroutine cont_eps_ab_ij(
     &              dist,
     &              irp_ab,vrt_a,eval_a,add_eval_a,
     &                     vrt_b,eval_b,add_eval_b,dim_dist,
     &              irp_ij,occ_i,eval_i,add_eval_i,
     &                     occ_j,eval_j,add_eval_j,
     &              energy
     &           )
      implicit none

c ARGUMENT LIST
      double precision dist(*)
      integer irp_ab, vrt_a(*), vrt_b(*)
      double precision eval_a(*), eval_b(*)
      integer add_eval_a(*), add_eval_b(*), dim_dist
      integer irp_ij, occ_i(*), occ_j(*)
      double precision eval_i(*), eval_j(*)
      integer add_eval_i(*), add_eval_j(*)
      double precision energy

c INTERNAL VARIABLES
      integer i, irp_i, max_i
      integer j, irp_j, max_j
      logical do_irp_j
      integer ndx
      double precision eps_j
      double precision eps_ij

c COMMON BLOCKS
#include "syminf.com" /* for nirrep and dirprd */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      ndx = 0
c   o assert irp_ij is in [1,nirrep]
      if ((irp_ij.lt.1).or.(nirrep.lt.irp_ij)) then
         print *, '@CONT_EPS_AB_IJ: Assertion failed.'
         print *, '   irp_ij = ',irp_ij
         ndx = -1
      end if
c   o assert irp_ab is in [1,nirrep]
      if ((irp_ab.lt.1).or.(nirrep.lt.irp_ab)) then
         print *, '@CONT_EPS_AB_IJ: Assertion failed.'
         print *, '   irp_ab = ',irp_ab
         ndx = -1
      end if
      if (ndx.ne.0) then
         call xcc_errex
      end if
#endif /* _ASSERT */

      ndx = 1

c   o loop over irreps of occupieds
      do irp_j = 1, nirrep
         irp_i = dirprd(irp_j,irp_ij)
         max_i = occ_i(irp_i)
         max_j = occ_j(irp_j)
         do_irp_j = ((max_j.ne.0).and.(max_i.ne.0))
      if (do_irp_j) then

c      o loop over orbitals IJ
         do j = 0, max_j-1
            eps_j  =         eval_j(add_eval_j(irp_j) + j)
         do i = 0, max_i-1
            eps_ij = eps_j + eval_i(add_eval_i(irp_i) + i)

c ----------------------------------------------------------------------

      call cont_eps_ab(
     &        dist(ndx),
     &        irp_ab,vrt_a,eval_a,add_eval_a,
     &               vrt_b,eval_b,add_eval_b,
     &        eps_ij,
     &        energy
     &     )

c ----------------------------------------------------------------------

c      o increment and close the loop over orbitals IJ
            ndx = ndx + dim_dist
c        end do i = 0, max_i-1
         end do
c        end do j = 0, max_j-1
         end do

c   o close the loop over irreps of occupieds
c     end if (do_irp_j)
      end if
c     end do irp_j = 1, nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

