
c This routine adds the quantity <0|adj(T2AA)*T1A^2|0> and its complex
c conjugate to the overlap. It can also be used for the T2BB piece if the
c argument ispin is 2; otherwise, if iuhf is 0, it doubles the T2AA term.

c INPUT
c    double precision dcore(*) : core space for doubles
c    double precision dscr(1)  : first free scratch address
c    integer          ispin    : Alpha/Beta signal (1=Alpha,2=Beta)

#include "disttype.h"

      subroutine ovrlp_t1t1t2aa(dcore,dscr,ispin)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), dscr(*)
      integer ispin

c INTERNAL VARIABLES
      integer irp_i, irp_j, irp_ij, irp_a, irp_b
      integer max_i, max_j, i, j, from_i, to_i, from_j, to_j, num_ij
      integer max_a, max_b, a, b, from_a, to_a, from_b, to_b, num_ab
      logical do_irp_j, do_irp_b, do_ovrlp, irp_i_eq_irp_j
      logical irp_a_eq_irp_b
      double precision t1_thresh, dtmp, ovrlp_tmp
      integer t2_start, t1_start, ndx, icount, itmp

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c The plan is:
c   o create a scratch vector with the same addressing as T2(IJ,AB)
c     for one whole I<J distribution (that's all A<B) made up of
c     T1(I,A)*T1(J,B)
c     NOTE: There exist T2s in which the symmetries of I and J
c           are not the same as those of A and B. The easiest example
c           is illustrated:
c              occ :  3-1-1-0
c              vrt : 10-3-6-2
c           The amplitude of A1,A1 -> B2,B2 is non-zero, but there is
c           no corresponding T1 amplitude A1 -> B2. Make sure to skip
c           these cases.
c   o ddot dscr with the T2 distribution and move on
c   o if this is an RHF calculation, double the contribution

      t2_start = ndx_t2(1,ispin)

c   o loop over the T2AA occupied irreps
      do irp_ij = 1,nirrep
      do irp_j  = 1,nirrep
         max_j  = pop(irp_j,ispin)
      if (max_j.ne.0) then
         irp_i  = dirprd(irp_j,irp_ij)
         max_i  = pop(irp_i,ispin)
         do_irp_j       = .false.
         irp_i_eq_irp_j = .false.
         if (irp_i.eq.irp_j) then
            irp_i_eq_irp_j = .true.
            if (max_j.gt.1) do_irp_j = .true.
         end if
         if (irp_i.lt.irp_j) then
            if (max_i.ne.0) do_irp_j = .true.
         end if
      if (do_irp_j) then

c      o loop over I<J occupied orbitals (as offsets)
         to_i   = max_i-1
         from_i = 0
         to_j   = max_j-1
         from_j = 0
         if (irp_i_eq_irp_j) from_j = 1
         do j = from_j,to_j
            if (irp_i_eq_irp_j) to_i = (j-1)
         do i = from_i,to_i

c         o loop over the T2AA virtual irreps
            do irp_b = 1,nirrep
               max_b = vrt(irp_b,ispin)
            if (max_b.ne.0) then
               irp_a = dirprd(irp_b,irp_ij)
               max_a = vrt(irp_a,ispin)
               do_irp_b = .false.
               irp_a_eq_irp_b = .false.
               if (irp_a.eq.irp_b) then
                  irp_a_eq_irp_b = .true.
                  if (max_b.gt.1) do_irp_b = .true.
               end if
               if (irp_a.lt.irp_b) then
                  if (max_a.ne.0) do_irp_b = .true.
               end if
            if (do_irp_b) then

c            o count the A<B pairs
               if (irp_a_eq_irp_b) then
c                  num_ab = rshift(max_b*(max_b-1),1)
                  num_ab = max_b*(max_b-1)/2
               else
                  num_ab = max_a*max_b
               end if

c            o do work only if the appropriate T1s exist
               do_ovrlp = .true.
c               if (irp_i.ne.irp_a) do_ovrlp = .false.
c               if (irp_j.ne.irp_b) do_ovrlp = .false.
               if (do_ovrlp) then

c ----------------------------------------------------------------------

c   o boundary offsets
      to_a   = max_a-1
      from_a = 0
      to_b   = max_b-1
      from_b = 0
      if (irp_a_eq_irp_b) from_b = 1

c   o debug this stupid thing
      itmp = 0
      do b = from_b,to_b
         if (irp_a_eq_irp_b) to_a = (b-1)
         do a = from_a,to_a
         write(*,'(1F17.12,3x,4(I1,A1,I2,A3))')
     &      dcore(t2_start+itmp),
     &      irp_a,'(',a+1,'), ',
     &      irp_b,'(',b+1,'), ',
     &      irp_i,'(',i+1,'), ',
     &      irp_j,'(',j+1,')  '
            itmp = itmp + 1
         end do
      end do

c ----------------------------------------------------------------------

c            o close the overlap flow-control loop
c              end if (do_ovrlp)
               end if

c            o get the address of the next A<B sublist
               t2_start = t2_start + num_ab

c         o close the T2AA virtual irrep flow-control loops
c           end if (do_irp_b)
            end if
c           end if (max_b.ne.0)
            end if
c           end do irp_b = 1,nirrep
            end do

c      o close the i<j orbital flow-control loops
c        end do i = from_i,to_i
         end do
c        end do j = from_j,to_j
         end do

c   o close the T2AA occupied irrep flow-control loops
c     end if (do_irp_j)
      end if
c     end if (max_j.ne.0)
      end if
c     end do irp_j  = 1,nirrep
      end do
c     end do irp_ij = 1,nirrep
      end do

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

