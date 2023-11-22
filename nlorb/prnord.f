      subroutine prnord(ordvc,nbas)
      implicit none
 
      integer nbas, ordvc(nbas), iii

      do iii = 1, nbas
C         write(*,10) iii, ordvc(iii)
      end do

 10   format('Molecular orbital ',I6,' became number ',I6)

      return
      end

