
#include "disttype.h"

      subroutine count(dcore)
      implicit none

c ARGUMENT LIST
      double precision dcore(*)

c INTERNAL VARIABLES
      integer irp_i, irp_j, irp_ij
      integer irp_a, irp_b
      integer i, j, max_i, max_j
      integer a, b, max_a, max_b, num_ab
      logical do_ovrlp
      double precision t1_thresh, dtmp, ovrlp_tmp
      integer t2_start, t1a_start, t1b_start, ndx, itmp, icount

c COMMON BLOCKS
#include "xcc.com"
#include "method.com"
#include "sym.com"

c ----------------------------------------------------------------------

#ifdef _PRINT_CRAP
      do i = 1, nirrep
         write(*,*) i,' : ',irpdpd(i,h_DistType_Ab),
     &                ' * ',irpdpd(i,h_DistType_Ij)
      end do
      write(*,*) 'ndx_t2(*,3) : ', (ndx_t2(i,3),i=1,nirrep)
      write(*,*) 'j : ',(pop(i,2),i=1,nirrep)
      write(*,*) 'i : ',(pop(i,1),i=1,nirrep)
      write(*,*) 'b : ',(vrt(i,2),i=1,nirrep)
      write(*,*) 'a : ',(vrt(i,1),i=1,nirrep)
#endif /* _PRINT_CRAP */

#ifdef _PRINT_T1_CRAP
      do irp_ij = 1, nirrep
         max_i = pop(irp_ij,1)
         max_a = vrt(irp_ij,1)
         do i = 0, max_i-1
         do a = 0, max_a-1
            write(*,'(1F16.10,3x,2(I1,A1,I2,A3))')
     &         dcore(ndx_t1(irp_ij,1)+i*max_a+a),
     &         irp_ij,'(',a+1,'), ',
     &         irp_ij,'(',i+1,')'
         end do
         end do
      end do
#endif /* _PRINT_T1_CRAP */

      icount = 0
      t2_start = ndx_t2(1,3)

c   o loop over the T2AB occupied irreps
      do irp_ij = 1, nirrep
      do irp_j  = 1, nirrep
         max_j  = pop(irp_j,2)
      if (max_j.ne.0) then
         irp_i  = dirprd(irp_j,irp_ij)
         max_i  = pop(irp_i,1)
      if (max_i.ne.0) then

         do j = 0, max_j-1
         do i = 0, max_i-1

      do irp_b = 1, nirrep
         max_b = vrt(irp_b,2)
      if (max_b.ne.0) then
         irp_a = dirprd(irp_b,irp_ij)
         max_a = vrt(irp_a,1)
      if (max_a.ne.0) then

c         do_ovrlp = .true.
c         if (irp_a.ne.irp_i) do_ovrlp = .false.
c         if (irp_b.ne.irp_j) do_ovrlp = .false.
c         if (do_ovrlp) then

         itmp = 0
         do b = 0, max_b-1
         do a = 0, max_a-1

c            write(*,*) a+1,b+1,i+1,j+1

            write(*,'(1F16.10,3x,4(I1,A1,I2,A3))')
     &         dcore(t2_start+itmp+a),
     &         irp_a,'(',a+1,'), ',
     &         irp_b,'(',b+1,'), ',
     &         irp_i,'(',i+1,'), ',
     &         irp_j,'(',j+1,')  '

            icount = icount + 1

         end do
         itmp = itmp + max_a
         end do

c        end if (do_ovrlp)
c         end if

         t2_start=t2_start+(max_a*max_b)

c   o close the irrep flow-control loops
c     end if (max_a.ne.0)
      end if
c     end if (max_b.ne.0)
      end if
c     end do irp_b = 1, nirrep
      end do

         end do
         end do

c     end if (max_i.ne.0)
      end if
c     end if (max_j.ne.0)
      end if
c     end do irp_j  = 1, nirrep
      end do
c     end do irp_ij = 1, nirrep
      end do

      write(*,*) 'There are ',icount,' T2AB amplitudes.'

c ----------------------------------------------------------------------

      return
      end

