      subroutine izeroarr(matrix,dim)
      implicit none  

      integer dim, matrix(dim,dim), iii, jjj

      do iii = 1, dim
         do jjj = 1, dim
            matrix(iii,jjj) = 0
         end do
      end do

      return
      end

