      subroutine crace(onehao,nbas)
      implicit none

      integer ibuf(600), int, luint, ilnbuf, nut,
     & indi, indj, nbas

      double precision onehao(nbas,nbas), buf(600)

      parameter(luint = 10)

      ilnbuf = 600

      open(luint,file = 'IIII    ',form = 'unformatted',
     & access = 'sequential')
      rewind luint     
 
      nut = ilnbuf 
      call locate(luint,'ONEHAMIL')
      do while (nut.eq.ilnbuf)
         read(luint) buf, ibuf, nut 
         do int = 1, nut
            call unpack2(ibuf(int),indi,indj)
            onehao(indi,indj) = buf(int)
            onehao(indj,indi) = buf(int)
         end do
      end do
    
      close(luint,status = 'KEEP')

      return
      end

