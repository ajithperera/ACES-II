
c This routine returns the quantity <0|adj(T2AA)*T1A^2|0>/2 and its complex
c conjugate for both spins. The functional form of this expression,
c i.e., the actual core of the subroutine, calculates:
c         n  [  N  (  AB      A B    A B   ) ]
c    2 * SUM [ SUM ( t   * [ t t  - t t  ] ) ]
c        I<J [ A<B (  IJ      I J    J I   ) ]

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision dscr(1)  : first free scratch address
c    double precision ovrlp    : the overlap

c OUTPUT
c    double precision ovrlp : the overlap plus the contribution
c                             from T2*T1T1/2

c#define _DEBUG_OVRLP_T2AA_T1T1

      subroutine ovrlp_t2aa_t1t1(dcore,dscr,ovrlp)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), dscr(*), ovrlp

c INTERNAL VARIABLES
      integer irp_i, irp_j, irp_ij
      integer irp_a, irp_b
      integer i, j, to_i, from_j,       max_i, max_j
      integer a, b, to_a, from_b, to_b, max_a, max_b, num_ab
      logical do_irp_j, do_irp_b
      double precision dtmp, ovrlp_tmp
      integer t2_start, t1_start, ndx, itmp, ispin

c PARAMETERS
c      double precision t1_thresh

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c   o loop over spin
      do ispin = 1, (1+iuhf)

c ----------------------------------------------------------------------

c The plan is:
c   o create a scratch vector with the same addressing as T2(IJ,AB)
c     for one whole I<J distribution (that's all A<B) made up of
c     T1(I,A)*T1(J,B)
c     NOTE: There exist T2 amplitudes in which the symmetries of I and J
c           are not the same as those of A and B. The easiest example
c           is illustrated:
c              occ :  3-1-1-0
c              vrt : 10-3-6-2
c           The amplitude of A1,A1 -> B2,B2 is non-zero, but there is
c           no corresponding T1 amplitude A1 -> B2. Make sure to skip
c           these cases.
c   o ddot dscr with the T2 distribution and repeat this for
c     T1(J,A)*T1(I,B) but subtract the dot product
c   o repeat for the next I<J pair
c   o double the contribution to account for the complex conjugate

      ovrlp_tmp = 0.0d0
      t2_start  = ndx_t2(1,ispin)

c   o loop over the T2AA occupied irreps
      do irp_ij = 1, nirrep
      do irp_j  = 1, nirrep
         max_j  = pop(irp_j,ispin)
      if (max_j.ne.0) then
         irp_i  = dirprd(irp_j,irp_ij)
         do_irp_j = .false.
         if (irp_i.lt.irp_j) then
            max_i = pop(irp_i,ispin)
            do_irp_j = (max_i.ne.0)
         end if
         if (irp_i.eq.irp_j) then
            max_i = max_j
            do_irp_j = (max_i.ne.1)
         end if
      if (do_irp_j) then

c      o loop over the I<J occupied orbitals (as offsets)
            to_i = max_i-1
            from_j = 0
            if (irp_i.eq.irp_j) from_j = 1
         do j = from_j, max_j-1
            if (irp_i.eq.irp_j) to_i = j-1
         do i = 0, to_i

c         o loop over the T2AA virtual irreps
            do irp_b = 1, nirrep
               max_b = vrt(irp_b,ispin)
            if (max_b.ne.0) then
               irp_a = dirprd(irp_b,irp_ij)
               do_irp_b = .false.
               if (irp_a.lt.irp_b) then
                  max_a = vrt(irp_a,ispin)
                  do_irp_b = (max_a.ne.0)
               end if
               if (irp_i.eq.irp_j) then
                  max_a = max_b
                  do_irp_b = (max_a.ne.1)
               end if
            if (do_irp_b) then

c            o count the A<B pairs
               if (irp_i.eq.irp_j) then
                  num_ab = rshift(max_b*(max_b-1),1)
               else
                  num_ab = max_a*max_b
               end if

c            o boundary offsets
               to_a = max_a-1
               to_b = max_b-1
               from_b = 0
               if (irp_i.eq.irp_j) from_b = 1

c ----------------------------------------------------------------------

c T1(I,A)*T1(J,B)

c   o do work only if the appropriate T1 amplitudes exist
      if ((irp_i.eq.irp_a).and.(irp_j.eq.irp_b)) then

c      o load T1(I,A) into the scratch vector
         ndx = 1
         t1_start = ndx_t1(irp_i,ispin) + max_a*i
         do b = from_b, to_b
            if (irp_i.eq.irp_j) to_a = b-1
            do a = 0, to_a
               dscr(ndx) = dcore(t1_start+a)
               ndx = ndx + 1
            end do
         end do

c      o scale the columns by T1(J,B)
         ndx = 1
         t1_start = ndx_t1(irp_j,ispin) + max_b*j
         do b = from_b, to_b
            dtmp = dcore(t1_start+b)
            if (irp_i.eq.irp_j) to_a = b-1
            do a = 0, to_a
               dscr(ndx) = dscr(ndx)*dtmp
               ndx = ndx + 1
            end do
         end do

c      o add T2(IJ,AB)*T1(I,A)*T1(J,B)
         do itmp = 0, num_ab-1
            ovrlp_tmp =   ovrlp_tmp
     &                  + (   dscr(1+itmp)
     &                      * dcore(t2_start+itmp) )
         end do

#ifdef _DEBUG_OVRLP_T2AA_T1T1
c      o debug
         itmp = 0
         do b = from_b, to_b
            dtmp = dcore(ndx_t1(irp_j,ispin) + max_b*j + b)
            if (irp_i.eq.irp_j) to_a = b-1
            do a = 0, to_a
               print '(3F16.10,3x,4(I1,A1,I2,A3))',
     &            dcore(t2_start+itmp), dtmp,
     &            dcore(ndx_t1(irp_i,ispin) + max_a*i + a),
     &            1+a,'[',irp_a,'], ',
     &            1+b,'[',irp_b,'], ',
     &            1+i,'[',irp_i,'], ',
     &            1+j,'[',irp_j,']  '
               itmp = itmp + 1
            end do
         end do
#endif /* _DEBUG_OVRLP_T2AA_T1T1 */

c   o close the overlap flow-control loop
c     end if ((irp_i.eq.irp_a).and.(irp_j.eq.irp_b))
      end if

c ----------------------------------------------------------------------

c T1(J,A)*T1(I,B)

c   o do work only if the appropriate T1 amplitudes exist
      if ((irp_j.eq.irp_a).and.(irp_i.eq.irp_b)) then

c      o load T1(J,A) into the scratch vector
         ndx = 1
         t1_start = ndx_t1(irp_j,ispin) + max_a*j
         do b = from_b, to_b
            if (irp_i.eq.irp_j) to_a = b-1
            do a = 0, to_a
               dscr(ndx) = dcore(t1_start+a)
               ndx = ndx + 1
            end do
         end do

c      o scale the columns by T1(I,B)
         ndx = 1
         t1_start = ndx_t1(irp_i,ispin) + max_b*i
         do b = from_b, to_b
            dtmp = dcore(t1_start+b)
            if (irp_i.eq.irp_j) to_a = b-1
            do a = 0, to_a
               dscr(ndx) = dscr(ndx)*dtmp
               ndx = ndx + 1
            end do
         end do

c      o negate T1(J,A)*T1(I,B)
         do itmp = 0, num_ab-1
            dscr(1+itmp) = -dscr(1+itmp)
         end do

c      o add T2(IJ,AB)*(-1)*T1(J,A)*T1(I,B)
         do itmp = 0, num_ab-1
            ovrlp_tmp =   ovrlp_tmp
     &                  + (   dscr(1+itmp)
     &                      * dcore(t2_start+itmp) )
         end do

#ifdef _DEBUG_OVRLP_T2AA_T1T1
c      o debug
         itmp = 0
         do b = from_b, to_b
            dtmp = dcore(ndx_t1(irp_i,ispin) + max_b*i + b)
            if (irp_i.eq.irp_j) to_a = b-1
            do a = 0, to_a
               print '(A,F15.10,2F16.10,3x,4(I1,A1,I2,A3))',
     &            '-',dcore(t2_start+itmp), dtmp,
     &            dcore(ndx_t1(irp_j,ispin) + max_a*j + a),
     &            1+a,'[',irp_a,'], ',
     &            1+b,'[',irp_b,'], ',
     &            1+i,'[',irp_i,'], ',
     &            1+j,'[',irp_j,']  '
               itmp = itmp + 1
            end do
         end do
#endif /* _DEBUG_OVRLP_T2AA_T1T1 */

c   o close the overlap flow-control loop
c     end if ((irp_j.eq.irp_a).and.(irp_i.eq.irp_b))
      end if

c ----------------------------------------------------------------------

c            o get the address of the next A<B T2(IJ,*) sublist
               t2_start = t2_start + num_ab

c         o close the T2AA virtual irrep flow-control loops
c           end if (do_irp_b)
            end if
c           end if (max_b.ne.0)
            end if
c           end do irp_b = 1, nirrep
            end do

c      o close the I<J orbital flow-control loops
c        end do i = from_i, to_i
         end do
c        end do j = from_j, to_j
         end do

c   o close the T2AA occupied irrep flow-control loops
c     end if (do_irp_j)
      end if
c     end if (max_j.ne.0)
      end if
c     end do irp_j  = 1, nirrep
      end do
c     end do irp_ij = 1, nirrep
      end do

c ----------------------------------------------------------------------

c   o add the complex conjugate
      ovrlp_tmp = ovrlp_tmp*2.0d0

      if (iuhf.eq.0) then
         ovrlp = ovrlp + ( 2.0d0 * ovrlp_tmp )
#ifdef _DEBUG_XCC
         print *, '<0|adj(T2AA)*T1A*T1A|0>/2 + c.c. = ',ovrlp_tmp
         print *, '<0|adj(T2BB)*T1B*T1B|0>/2 + c.c. = ',ovrlp_tmp
#endif /* _DEBUG_XCC */
      else
         ovrlp = ovrlp + ovrlp_tmp
#ifdef _DEBUG_XCC
         if (ispin.eq.1) then
            print *, '<0|adj(T2AA)*T1A*T1A|0>/2 + c.c. = ',ovrlp_tmp
         else
            print *, '<0|adj(T2BB)*T1B*T1B|0>/2 + c.c. = ',ovrlp_tmp
         end if
#endif /* _DEBUG_XCC */
      end if

c ----------------------------------------------------------------------

c     end do ispin = 1, (1+iuhf)
      end do

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

