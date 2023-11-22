      subroutine dzeroarr(dim1,dim2,arr)
      implicit none

      integer dim1, dim2, iii, jjj
 
      double precision arr(dim1,dim2)

      do iii = 1, dim1
         do jjj = 1, dim2
            arr(iii,jjj) = 0.D0
         end do
      end do

      return
      end
     
