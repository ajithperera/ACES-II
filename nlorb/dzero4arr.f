      subroutine dzero4arr(dim1,dim2,dim3,dim4,arr)
      implicit none

      integer dim1, dim2, dim3, dim4, iii, jjj, kkk, lll
 
      double precision arr(dim1,dim2,dim3,dim4)

      do iii = 1, dim1
         do jjj = 1, dim2
            do kkk = 1, dim3
               do lll = 1, dim4
                  arr(iii,jjj,kkk,lll) = 0.D0
               end do
            end do
         end do
      end do

      return
      end
     
