
c This routine prints a header full of information.

      subroutine print_header
      implicit none

c INTERNAL VARIABLES
      integer i

c COMMON BLOCKS
#include "sym.com" /* for pop(,) and vrt(,) */
#include "syminf.com" /* for nirrep */

c ----------------------------------------------------------------------

      print '(/)'
      print '(24x,a)', '********************************'
      print '(24x,a)', '* THE FASTEST TAU2 IN THE WEST *'
      print '(24x,a)', '********************************'
      print '(/)'

      print '(22x,a)', '         --OCCUPIED--   --VIRTUALS--'
      print '(22x,a)', 'irrep    alpha   beta   alpha   beta'
      print '(22x,a)', '-----    -----   ----   -----   ----'
      do i = 1, nirrep
         print    '(22x,i5,  4x, i5,  3x,i4, 3x,i5,  3x,i4)',
     &                  i, pop(i,1), pop(i,2), vrt(i,1), vrt(i,2)
      end do
      print '(/)'

      return
c     end subroutine print_header
      end

