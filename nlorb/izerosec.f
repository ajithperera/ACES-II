      subroutine izerosec(nchoose,bondsize)
      implicit none
  
      integer bondsize, nchoose(bondsize), iii

      do iii = 1, bondsize
         nchoose(iii) = 0
      end do

      return
      end 

