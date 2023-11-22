      subroutine scalmat(matrix,dim,num)
      implicit none

      integer dim, iii, jjj 

      double precision num, matrix(dim,dim)

      do iii = 1, dim
         do jjj = 1, dim
            matrix(iii,jjj) = num*matrix(iii,jjj)
         end do
      end do
 
      return
      end

