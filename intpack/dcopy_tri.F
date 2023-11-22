      subroutine dcopy_tri(lenght,old,new)
#include <aces.h>
       
      integer lenght

      M_REAL old (lenght,lenght,lenght),
     &       new (lenght,lenght,lenght)

      integer i,j,k
      
      do i=1,lenght
         do j=1,lenght
            do k=1,lenght
               new(i,j,k)=old(i,j,k)
            end do
         end do
      end do
      return
      end
