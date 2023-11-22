
c This routine returns the quantity <0|adj(T2AB)*T1A*T1B|0>/2 and its
c complex conjugate.

c INPUT
c    double precision dcore(*) : address space for xcc.com
c    double precision dscr(*)  : first free scratch address
c    double precision ovrlp    : the overlap

c OUTPUT
c    double precision ovrlp : the overlap plus the contribution
c                             from T2*T1T1/2

c#define _DEBUG_OVRLP_T2AB_T1T1

      subroutine ovrlp_t2ab_t1t1(dcore,dscr,ovrlp)
      implicit none

c ARGUMENT LIST
      double precision dcore(*), dscr(*), ovrlp

c INTERNAL VARIABLES
      integer irp_i, irp_j, irp_ij
      integer irp_a, irp_b
      integer i, j, max_i, max_j
      integer a, b, max_a, max_b
      double precision dtmp, ovrlp_tmp
      integer t2_start, t1a_start, t1b_start, ioff

c PARAMETERS
      double precision zero, two, t1_thresh
      parameter (zero=0.0d0, two=2.0d0, t1_thresh=1.0d-20)

c COMMON BLOCKS
#include "xcc.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Entered file ',__FILE__
#endif

c ----------------------------------------------------------------------

c The plan is:
c   o loop over I, j, and b in T2(Ij,Ab)
c      - create a scratch vector of T1(I,A)*T1(j,b) for all A
c      - ddot dscr with the T2 distribution and move on

      ovrlp_tmp = zero
      t2_start  = ndx_t2(1,3)

c   o loop over occupied irreps and orbitals in T2(Ij,Ab)
      do irp_ij = 1, nirrep
      do irp_j  = 1, nirrep
         max_j  = pop(irp_j,2)
      if (max_j.ne.0) then
         irp_i  = dirprd(irp_j,irp_ij)
         max_i  = pop(irp_i,1)
      if (max_i.ne.0) then
      do j = 0, max_j-1
      do i = 0, max_i-1

c      o loop over virtual irreps in T2(Ij,Ab)
         do irp_b = 1, nirrep
            max_b = vrt(irp_b,2)
         if (max_b.ne.0) then
            irp_a = dirprd(irp_b,irp_ij)
            max_a = vrt(irp_a,1)
         if (max_a.ne.0) then

c         o filter-out improper symmetry
            if ((irp_a.eq.irp_i).and.(irp_b.eq.irp_j)) then

c ----------------------------------------------------------------------

c   o initialize the offsets and addresses
      ioff = 0
      t1b_start = ndx_t1(irp_j,2) + max_b*j
      t1a_start = ndx_t1(irp_i,1) + max_a*i

c   o dot T1(I,A) with T2(Ij,Ab) and scale by T1(j,b)
      do b = 0, max_b-1
         if (t1_thresh.lt.abs(dcore(t1b_start+b))) then
            dtmp = zero
            do a = 0, max_a-1
               dtmp =   dtmp
     &                + (   dcore(t1a_start+a)
     &                    * dcore(t2_start+ioff+a) )
            end do
            ovrlp_tmp = ovrlp_tmp + ( dtmp * dcore(t1b_start+b) )
#ifdef _DEBUG_OVRLP_T2AB_T1T1
            do a = 0, max_a-1
               print '(3F16.10,3x,4(I1,A1,I2,A3))',
     &            dcore(t2_start+ioff+a),
     &            dcore(t1b_start+b),
     &            dcore(t1a_start+a),
     &            1+a,'[',irp_a,'], ',
     &            1+b,'[',irp_b,'], ',
     &            1+i,'[',irp_i,'], ',
     &            1+j,'[',irp_j,']  '
            end do
#endif /* _DEBUG_OVRLP_T2AB_T1T1 */
c        end if (t1_thresh.lt.abs(dcore(t1b_start+b)))
         end if
         ioff = ioff + max_a
      end do

c ----------------------------------------------------------------------

c         o close the symmetry conditional
c           end if ((irp_a.eq.irp_i).and.(irp_b.eq.irp_j))
            end if

c         o point to the next block of virtual irreps
            t2_start = t2_start + (max_a*max_b)

c      o close the loop over virtual irreps in T2(Ij,Ab)
c        end if (max_a.ne.0)
         end if
c        end if (max_b.ne.0)
         end if
c        end do irp_b = 1, nirrep
         end do

c   o close the loop over occupied irreps and orbitals in T2(Ij,Ab)
c     end do i = 0, max_i-1
      end do
c     end do j = 0, max_j-1
      end do
c     end if (max_i.ne.0)
      end if
c     end if (max_j.ne.0)
      end if
c     end do irp_j  = 1, nirrep
      end do
c     end do irp_ij = 1, nirrep
      end do

c ----------------------------------------------------------------------

c   o add the complex conjugate
      ovrlp = ovrlp + ( two * ovrlp_tmp )

#ifdef _DEBUG_XCC
      print *, '<0|adj(T2AB)*T1A*T1B|0>/2 + c.c. = ',ovrlp_tmp*two
#endif /* _DEBUG_XCC */

c ----------------------------------------------------------------------

#ifdef _TRACE_XCC
      print *, '@XCC: Leaving file ',__FILE__
#endif

      return
      end

