      subroutine cartdeg(degcart,angmax)
      implicit none

      integer angmax, degcart(angmax), iii,
     & temp

      temp = 0

      do iii = 1, angmax
         degcart(iii) = temp + iii
         temp = degcart(iii)
      end do 

      return
      end

