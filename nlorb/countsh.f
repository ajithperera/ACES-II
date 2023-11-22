      subroutine countsh(natoms,angmax,count,count2,degen)
      implicit none

      integer natoms, angmax, count(natoms,angmax),
     & count2(natoms,angmax), degen(angmax), iii, lll

      do iii = 1, natoms
         do lll = 1, angmax
            count2(iii,lll) = count(iii,lll)/degen(lll)
         end do
      end do
 
      return
      end

