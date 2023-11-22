
c This routine adds the quantity <0|adj(T2AB)*T1B^2|0> and its complex
c conjugate to the overlap.

c INPUT
c    double precision dcore(*) : core space for doubles
c    double precision dscr(1)  : first free scratch address

#include "disttype.h"

      subroutine ovrlp_t1t1t2ab(dcore,dscr)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), dscr(*)

c INTERNAL VARIABLES
      integer irp_i, irp_j, irp_ij, irp_a, irp_b
      integer max_i, max_j, i, j, from_i, to_i, from_j, to_j, num_ij
      integer max_a, max_b, a, b, from_a, to_a, from_b, to_b, num_ab
      logical do_irp_j, do_irp_b, do_ovrlp, irp_i_eq_irp_j
      double precision t1_thresh, dtmp, ovrlp_tmp
      integer t2_start, t1a_start, t1b_start, ndx, icount, itmp

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
c   o create a scratch vector with the same addressing as T2(Ij,Ab)
c     for one whole b;Ij distribution (that's all A) made up of
c     T1(I,A)*T1(j,b)
c   o ddot dscr with the T2 distribution and move on

      ovrlp_tmp = 0.0d0
      t2_start = ndx_t2(1,3)

c   o loop over the T2AB occupied irreps
      do irp_ij = 1,nirrep
      do irp_j  = 1,nirrep
         max_j  = pop(irp_j,2)
      if (max_j.ne.0) then
         irp_i  = dirprd(irp_j,irp_ij)
         max_i  = pop(irp_i,1)
         do_irp_j = .false.
         if (max_i.ne.0) do_irp_j = .true.
      if (do_irp_j) then

c      o loop over Ij occupied orbitals (as offsets)
         to_i   = max_i-1
         from_i = 0
         to_j   = max_j-1
         from_j = 0
         do j = from_j,to_j
         do i = from_i,to_i

c         o loop over the T2AA virtual irreps
            do irp_b = 1,nirrep
               max_b = vrt(irp_b,2)
            if (max_b.ne.0) then
               irp_a = dirprd(irp_b,irp_ij)
               max_a = vrt(irp_a,1)
               do_irp_b = .false.
               if (max_a.ne.0) do_irp_b = .true.
            if (do_irp_b) then

c            o count the Ab pairs
               num_ab = max_a*max_b

c            o do work only if the appropriate T1s exist
               do_ovrlp = .true.
               if (irp_i.ne.irp_a) do_ovrlp = .false.
               if (irp_j.ne.irp_b) do_ovrlp = .false.
               if (do_ovrlp) then

c ----------------------------------------------------------------------

c   o boundary offsets
      to_a   = max_a-1
      from_a = 0
      to_b   = max_b-1
      from_b = 0

c   o for each b, load T1(I,A), scale by T1(j,b), and dot with T2(Ij,Ab)
      t1a_start = ndx_t1(irp_i,1) + max_a*i
      t1b_start = ndx_t1(irp_j,2) + max_b*j
      do b = from_b,to_b
         do a = from_a,to_a
            dscr(1+a) = dcore(t1a_start+a)
         end do
         dtmp = dcore(t1b_start+b)
         do a = from_a,to_a
            dscr(1+a) = dscr(1+a)*dtmp
         end do
         itmp = b*max_b
         do a = from_a,to_a
            ovrlp_tmp =   ovrlp_tmp
     &                  + (   dscr(1+a)
     &                      * dcore(t2_start+itmp+a) )
         end do
      end do

c   o debug this stupid thing
c      itmp = 0
c      do b = from_b,to_b
c         dtmp = dcore(ndx_t1(irp_j,ispin) + max_b*j + b)
c         if (irp_i_eq_irp_j) to_a = (b-1)
c         do a = from_a,to_a
c         write(*,'(2F17.12,3x,4(I1,A1,I2,A3))')
c     &      dcore(t2_start+itmp),
c     &      dtmp * dcore(ndx_t1(irp_i,ispin) + max_a*i + a),
c     &      irp_a,'(',a+1,'), ',
c     &      irp_b,'(',b+1,'), ',
c     &      irp_i,'(',i+1,'), ',
c     &      irp_j,'(',j+1,')  '
c            itmp = itmp + 1
c         end do
c      end do

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

      ovrlp_tmp = 2.0d0*ovrlp_tmp
      overlap = overlap + ovrlp_tmp
      write(*,*) '@XCC: <0|adj(T2AB)*T1A*T1B|0> = ',ovrlp_tmp

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      write(*,*) '@XCC: Leaving file ',__FILE__
#endif

      return
      end

