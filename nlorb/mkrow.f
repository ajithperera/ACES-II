      subroutine mkrow(row,numrow,angbf,angmax,nbas,nbastot)
      implicit none
 
      integer numrow, row(numrow), angbf(nbas), angmax, nbas,
     & nbastot, degcart(angmax), degsph(angmax), iii, count,
     & rownum

      call sphdeg(degsph,angmax)
      call cartdeg(degcart,angmax)
 
      count = 0
      rownum = 0
      do iii = 1, nbas
         if (angbf(iii).ge.2) then
            count = count + 1
         end if
         if (count.gt.degsph(angbf(iii)+1)) then
            rownum = rownum + 1
            row(rownum) = iii
         end if
         if (count.eq.degcart(angbf(iii)+1)) then
            count = 0
         end if
      end do

      return
      end

