
      subroutine iprt(nrow,ncol,mat,title)
#include <aces.h>
      integer nrow,ncol,mat(nrow,ncol)
      character *(*) title
      integer i,j
      write (6,*) '*** ',title
      do i=1,nrow
        write(6,100) (mat(i,j),j=1,ncol)
      end do
      write(6,*)
  100 format (15(i4,1x))
      return
      end
