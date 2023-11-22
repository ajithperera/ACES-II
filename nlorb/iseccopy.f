      subroutine iseccopy(inlst,outlst,dim)
      implicit none

      integer dim, iii, inlst(dim), outlst(dim)

      do iii = 1, dim
         outlst(iii) = inlst(iii)
      end do
 
      return
      end
 
