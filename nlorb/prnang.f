      subroutine prnang(largel,nbas,mombf,centbf)
      implicit none

      integer iii

      integer largel, nbas

      integer mombf(nbas), centbf(nbas)

      largel = 0

      do iii = 1, nbas
         write(*,10) iii, centbf(iii), mombf(iii)
         if (largel.lt.mombf(iii)) then
            largel = mombf(iii)
         end if
      end do 

      write(*,*)
      write(*,*)

 10   format ('Function ',I5,' on center ',I5,' with ang mom ',I5)

      return
      end

