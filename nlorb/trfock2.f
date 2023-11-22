      subroutine trfock2(onehao,geemx,nbas,fock)
      implicit none

      integer iii, jjj, nbas
      double precision onehao(nbas,nbas), 
     & geemx(nbas,nbas), fock(nbas,nbas)

      do iii = 1, nbas
         do jjj = 1, nbas 
            fock(iii,jjj) = onehao(iii,jjj) + geemx(iii,jjj)
         end do
      end do

      return
      end

