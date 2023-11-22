      subroutine unpkoneh(onehpk,onehunpk,nbas)
      implicit none

      integer muu, nhu, index

      integer nbas

      double precision onehpk(nbas*(nbas+1)/2),
     & onehunpk(nbas,nbas)

      do index = 1, nbas*(nbas+1)/2
         call unpack2(index,muu,nhu) 
         onehunpk(muu,nhu) = onehpk(index)
         onehunpk(nhu,muu) = onehpk(index)
      end do 

      return
      end

