      subroutine dtranspose(matin,matout,dim)
      implicit none
 
      integer dim, iii, jjj
      double precision matin(dim,dim), matout(dim,dim)

      do iii = 1, dim
         do jjj = 1, dim
            matout(iii,jjj) = matin(jjj,iii)
         end do
      end do
 
      return
      end

