      subroutine countmom(natoms,nbas,centbf,mombf,angmax,count)
      implicit none

      integer natoms, nbas, centbf(nbas), mombf(nbas), angmax,
     & count(natoms,angmax), iii, jjj, lll

      do iii = 1, natoms
         do jjj = 1, nbas
            if (centbf(jjj).eq.iii) then
               do lll = 1, angmax
                  if (mombf(jjj).eq.(lll - 1)) then
                     count(iii,lll) = count(iii,lll) + 1
                  end if
               end do
            end if
         end do
      end do
  
      return
      end


       
