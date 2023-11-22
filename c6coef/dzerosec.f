      subroutine dzerosec(len,vec)
      implicit none
  
      integer len, iii
  
      double precision vec(len)

      do iii = 1, len
         vec(iii) = 0.D0
      end do

      return
      end 

