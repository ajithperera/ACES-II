      subroutine delrow(row,numrow,wrow,worow,nbas)
      implicit none

      integer wrow(nbas,nbas), worow(nbas-numrow,nbas),
     & numrow, nbas, row(numrow), count, iii, jjj, kkk,
     & badrow 

      count = 0
      badrow = 0

      do iii = 1, nbas

         do kkk = 1, numrow
            if (iii.eq.row(kkk)) then
               count = count + 1
               badrow = 1
            end if
         end do
       
         do jjj = 1, nbas
            if (badrow.ne.1) then
               worow(iii-count,jjj) = wrow(iii,jjj)
            end if
         end do

      badrow = 0
      end do
 
      return
      end
 
