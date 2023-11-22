      subroutine totalden(dena,denb,total,dim)
      implicit none

      integer dim, iii, jjj
    
      double precision dena(dim,dim), denb(dim,dim),
     & total(dim,dim)
 
      do iii = 1, dim
         do jjj = 1, dim
            total(iii,jjj) = dena(iii,jjj) + denb(iii,jjj)
         end do
      end do
 
      return
      end

