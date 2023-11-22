      subroutine orthock(ovrmo,ovrao,coef,nbas)
      implicit none

      integer nbas, iii, jjj, mu, nhu
      double precision ovrmo(nbas,nbas), ovrao(nbas,nbas),
     &                 coef(nbas,nbas)

      do iii = 1, nbas
         do jjj = 1, nbas
            do mu = 1, nbas
               do nhu = 1, nbas
                  ovrmo(iii,jjj) = ovrmo(iii,jjj) +
     & coef(mu,iii)*coef(nhu,jjj)*ovrao(mu,nhu)
               end do
	    end do
         end do
      end do

      return
      end

