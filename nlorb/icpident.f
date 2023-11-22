      subroutine icpident(matrix,dim)
      implicit none

      integer dim, matrix(dim,dim), iii
  
      do iii = 1, dim
         matrix(iii,iii) = 1
      end do
 
      return
      end

