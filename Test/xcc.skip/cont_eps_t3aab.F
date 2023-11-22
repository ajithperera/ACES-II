
c This routine divides a T3(IJk,ABc) distribution of pre-amplitudes by
c epsilon(IJk,ABc) thus creating pure T3 amplitudes. It also calculates
c the E[4]T contribution from the distribution and the T3 overlap.

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision t3aab(*) : the IJk distribution of T3 pre-amplitudes
c    double precision eps_ijk  : the first half of the epsilon denominator
c    integer          irp_ijk  : the total symmetry of the distribution
c    integer          ispin1   : the primary   spin case to handle
c    integer          ispin2   : the secondary spin case to handle
c    double precision e4t      : the E[4]T correlation energy (not zeroed)
c    double precision t3_ovrlp : <0|adj(T3AAB)*T3AAB|0> (not zeroed)

      subroutine cont_eps_t3aab(dcore,t3aab,
     &                          eps_ijk,irp_ijk,
     &                          ispin1,ispin2,
     &                          e4t,t3_ovrlp)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), t3aab(*), eps_ijk, e4t, t3_ovrlp
      integer irp_ijk, ispin1, ispin2

c INTERNAL VARIABLES
      integer a, irp_a, max_a,         to_a, start_a
      integer b, irp_b, max_b, from_b,       start_b
      integer c, irp_c, max_c,               start_c
      integer irp_bc
      logical do_irp_b
      integer ndx
      double precision eps_ijk_c
      double precision eps_ijk_bc
      double precision eps_ijk_abc
      double precision dtmp, tmp_e4t, tmp_ovrlp

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

#ifdef _ASSERT
      ndx = 0
c   o assert ispin1 and ispin2 are 1 or 2
      if ((ispin1.lt.1).or.(ispin1.gt.2).or.
     &    (ispin2.lt.1).or.(ispin2.gt.2)    ) then
         print *, '@CONT_EPS_T3AAB: Assertion failed.'
         print *, '   ispin1 = ',ispin1
         print *, '   ispin2 = ',ispin2
         ndx = -1
      end if
c   o assert irp_ijk is in [1,nirrep]
      if ((irp_ijk.lt.1).or.(nirrep.lt.irp_ijk)) then
         print *, '@CONT_EPS_T3AAB: Assertion failed.'
         print *, '   irp_ijk = ',irp_ijk
         ndx = -1
      end if
      if (ndx.ne.0) then
         call xcc_errex
      end if
#endif /* _ASSERT */

      ndx = 1
      tmp_ovrlp = 0.0d0
      tmp_e4t   = 0.0d0

c   o loop over irreps of virtuals in T3(IJk,ABc)
      do irp_c = 1, nirrep
         max_c = vrt(irp_c,ispin2)
      if (max_c.ne.0) then
      do irp_b = 1, nirrep
         irp_bc = dirprd(irp_b,irp_c)
         irp_a  = dirprd(irp_bc,irp_ijk)
         if (irp_a.lt.irp_b) then
            max_a = vrt(irp_a,ispin1)
            max_b = vrt(irp_b,ispin1)
            do_irp_b = ((max_b.ne.0).and.(max_a.ne.0))
         else
            if (irp_a.gt.irp_b) then
               do_irp_b = .false.
            else
               max_a = vrt(irp_a,ispin1)
               max_b = max_a
               do_irp_b = (max_b.gt.1)
            end if
         end if
      if (do_irp_b) then

c      o retrieve starting addresses once
         start_a = ndx_eval(irp_a,2,ispin1)
         start_b = ndx_eval(irp_b,2,ispin1)
         start_c = ndx_eval(irp_c,2,ispin2)

c      o set orbital boundary conditions (as offsets)
         to_a = max_a-1
         if (irp_a.lt.irp_b) then
            from_b = 0
         else
            from_b = 1
         end if

c      o loop over orbitals A<B,c
         do c = 0,      max_c-1
            eps_ijk_c   = eps_ijk    - dcore(start_c + c)
         do b = from_b, max_b-1
            eps_ijk_bc  = eps_ijk_c  - dcore(start_b + b)
            if (irp_a.eq.irp_b) to_a = b-1
         do a = 0,      to_a
            eps_ijk_abc = eps_ijk_bc - dcore(start_a + a)

c ----------------------------------------------------------------------

      t3aab(ndx) = t3aab(ndx) / eps_ijk_abc
      dtmp = t3aab(ndx) * t3aab(ndx)
      tmp_ovrlp = tmp_ovrlp +   dtmp
      tmp_e4t   = tmp_e4t   + ( dtmp * eps_ijk_abc )

c ----------------------------------------------------------------------

c         o increment the T3 pointer
            ndx = ndx + 1

c      o close the loop over orbitals A<B,c
c        end do a = 0, to_a
         end do
c        end do b = from_b, max_b-1
         end do
c        end do c = 0, max_c-1
         end do

c   o close the loop over irreps of virtuals in T3(IJk,ABc)
c     end if (do_irp_b)
      end if
c     end do irp_b = 1, nirrep
      end do
c     end if (max_c.ne.0)
      end if
c     end do irp_c = 1, nirrep
      end do

c   o add in the T3 overlap and E[4]T contributions
      t3_ovrlp = t3_ovrlp + tmp_ovrlp
      e4t      = e4t      + tmp_e4t

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

