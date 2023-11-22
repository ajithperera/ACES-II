      subroutine eneval2(denmxtot,dena,denb,onehao,
     & focka,fockb,nbas,eelec)
      implicit none

      integer nbas, iii, jjj

      double precision denmxtot(nbas,nbas), dena(nbas,nbas),
     & denb(nbas,nbas), onehao(nbas,nbas), focka(nbas,nbas),
     & fockb(nbas,nbas), eelec

      eelec = 0.0D0
 
      do iii = 1, nbas
         do jjj = 1, nbas
            eelec = eelec + 0.5D0*(denmxtot(jjj,iii)*onehao(iii,jjj) +
     & dena(jjj,iii)*focka(iii,jjj) + denb(jjj,iii)*fockb(iii,jjj))
         end do
      end do

      return
      end
