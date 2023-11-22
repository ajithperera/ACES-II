      subroutine dcopyarr(inarr,outarr,dim)
      implicit none
 
      integer dim, iii, jjj
 
      double precision inarr(dim,dim), outarr(dim,dim)

      do iii = 1, dim
         do jjj = 1, dim
            outarr(iii,jjj) = inarr(iii,jjj)
         end do
      end do

      return
      end

