      subroutine eneval(denmx,onehao,fock,nbas,eelec)
      implicit none

      integer nbas, iii, jjj
      double precision denmx(nbas,nbas),
     & onehao(nbas,nbas), fock(nbas,nbas), eelec

      eelec = 0.0D0
 
      do iii = 1, nbas
         do jjj = 1, nbas
            eelec = eelec + 0.5D0*denmx(jjj,iii)*(onehao(iii,jjj) +
     & fock(iii,jjj))
         end do
      end do

      return
      end
