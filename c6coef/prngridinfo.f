      subroutine prngridinfo(numpoints,wghts,freqs,grdvl)
      implicit none

      integer ind, numpoints

      double precision wghts(numpoints), freqs(numpoints),
     & grdvl(numpoints)

      do ind = 1, numpoints
         write(*,10) ind, grdvl(ind), freqs(ind), wghts(ind)
      end do

 10   format(I6,3F20.15)

      return
      end
