      subroutine mkdens(coef,denmx,ordvc,nbas,nocc,ordered)
      implicit none

      integer nbas, nelec, nocc, iii, jjj, kkk, lll,
     & ordvc(nbas)

      double precision coef(nbas,nbas), denmx(nbas,nbas)

      logical ordered

      do iii = 1, nbas
         do jjj = 1, nbas
            do kkk = 1, nocc 
               lll = kkk
               if (ordered) then
                  lll = ordvc(kkk)
               end if
               denmx(iii,jjj) = denmx(iii,jjj) +
     & coef(iii,lll)*coef(jjj,lll)
            end do
         end do
      end do

      return
      end

