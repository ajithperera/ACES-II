      subroutine doeval(fockmo,evals,nbastot)
      implicit none

      integer nbastot, iii
      double precision fockmo(nbastot,nbastot), evals(nbastot)

      do iii = 1, nbastot
         evals(iii) = fockmo(iii,iii)
      end do
 
      return
      end

