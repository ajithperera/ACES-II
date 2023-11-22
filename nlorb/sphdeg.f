      subroutine sphdeg(degen,angmax)
      implicit none
 
      integer angmax, degen(angmax), iii

      do iii = 1, angmax
         degen(iii) = 2*(iii - 1) + 1
      end do

      return
      end

