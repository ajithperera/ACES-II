
c This routine divides one operator symmetry block by epsilon. Specifically,
c the distribution must be stored all A<B where A and B span the same virtual
c space for all I<J where I and J span the same occupied space. Additionally,
c the quantity O**2/eps is returned. For example, if O is the 2-electron
c antisymmetrized integral, then the energy returned is a contribution to
c the MBPT(2) energy.

c INPUT:
c    double precision dist(*)  : the distributions to modify

c    integer irp_ab : the symmetry of the virtual distributions
c    integer vrt_a(*) : the population vector of the virtual irreps
c    double precision eval_a(*)     : the virtual eigenvalues
c    integer          add_eval_a(*) : the starting addresses by irrep
c                                     of virtual eigenvalues in eval_a(*)
c    integer dim_dist : the length of one A<B distribution

c    integer irp_ij : the symmetry of the occupied distributions
c    integer occ_i(*) : the population vector of the occupied irreps
c    double precision eval_i(*)     : the occupied eigenvalues
c    integer          add_eval_i(*) : the starting addresses by irrep
c                                     of occupied eigenvalues in eval_i(*)

c OUTPUT:
c    double precision energy : the O**2/eps contribution (NOT ZEROED)

      subroutine cont_eps_alta_ilti(
     &              dist,
     &              irp_ab,vrt_a,eval_a,add_eval_a,dim_dist,
     &              irp_ij,occ_i,eval_i,add_eval_i,
     &              energy
     &           )
      implicit none

c ARGUMENT LIST
      double precision dist(*)
      integer irp_ab, vrt_a(*)
      double precision eval_a(*)
      integer add_eval_a(*), dim_dist
      integer irp_ij, occ_i(*)
      double precision eval_i(*)
      integer add_eval_i(*)
      double precision energy

c INTERNAL VARIABLES
      integer i, irp_i, max_i,         to_i
      integer j, irp_j, max_j, from_j
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
         print *, '@CONT_EPS_ALTA_ILTI: Assertion failed.'
         print *, '   irp_ij = ',irp_ij
         ndx = -1
      end if
c   o assert irp_ab is in [1,nirrep]
      if ((irp_ab.lt.1).or.(nirrep.lt.irp_ab)) then
         print *, '@CONT_EPS_ALTA_ILTI: Assertion failed.'
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
         if (irp_i.lt.irp_j) then
            max_i = occ_i(irp_i)
            max_j = occ_i(irp_j)
            do_irp_j = ((max_j.ne.0).and.(max_i.ne.0))
         else
            if (irp_i.gt.irp_j) then
               do_irp_j = .false.
            else
               max_i = occ_i(irp_i)
               max_j = max_i
               do_irp_j = (max_j.gt.1)
            end if
         end if
      if (do_irp_j) then

c      o set orbital boundary conditions (as offsets)
         to_i = max_i-1
         if (irp_i.lt.irp_j) then
            from_j = 0
         else
            from_j = 1
         end if

c      o loop over orbitals I<J
         do j = from_j, max_j-1
            eps_j  =         eval_i(add_eval_i(irp_j) + j)
            if (irp_i.eq.irp_j) to_i = j-1
         do i = 0, to_i
            eps_ij = eps_j + eval_i(add_eval_i(irp_i) + i)

c ----------------------------------------------------------------------

      call cont_eps_alta(
     &        dist(ndx),
     &        irp_ab,vrt_a,eval_a,add_eval_a,
     &        eps_ij,
     &        energy
     &     )

c ----------------------------------------------------------------------

c      o increment and close the loop over orbitals I<J
            ndx = ndx + dim_dist
c        end do i = 0, to_i
         end do
c        end do j = from_j, max_j-1
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

