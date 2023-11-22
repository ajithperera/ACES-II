      subroutine mkgee4(geemx,nbas,gal,denmx)
      implicit none

      integer iii, jjj, kkk, lll, nbas
      double precision geemx(nbas,nbas), denmx(nbas,nbas),
     & gal(nbas,nbas,nbas,nbas) 

      do iii = 1, nbas
	 do jjj = 1, nbas
            do kkk = 1, nbas 
               do lll = 1, nbas 
                  geemx(iii,jjj) = geemx(iii,jjj) + denmx(kkk,lll)*
     & (gal(iii,jjj,kkk,lll)-0.5D0*gal(iii,lll,kkk,jjj))
               end do
            end do
	 end do
      end do

      return
      end

