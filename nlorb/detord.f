      subroutine detord(owntrans,ordvec,nbas)
      implicit none

      integer nbas, ordvec(nbas), iii, jjj
      
      double precision owntrans(nbas,nbas), thres

      parameter (thres = 0.D-12)

      do jjj = 1, nbas
         do iii = 1, nbas
            if (abs(owntrans(iii,jjj)-1.D0).le.thres) then
               if (iii.ne.jjj) then
C                  write(*,10) jjj, iii
                  ordvec(jjj) = iii
               end if
               if (iii.eq.jjj) then
C                  write(*,11) jjj, iii
                  ordvec(jjj) = iii
               end if
            end if
         end do
      end do
 
 10   format('Offdiag: old ',I6,' is new ',I6)
 11   format('Diag: old ',I6,' is new ',I6)

      return
      end
