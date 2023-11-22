      subroutine ovraces(ovrlp,nbas)
      implicit none

      integer ibuf(600), ilnbuf, luint, nut, int,
     & indi, indj, nbas

      double precision ovrlp(nbas,nbas), buf(600)

      parameter(luint = 10)

      ilnbuf = 600

      open(luint,file = 'IIII    ',form = 'unformatted',
     & access = 'sequential') 
      rewind luint
      
      nut = ilnbuf
      call locate(luint,'OVERLAP ')
      do while (nut.eq.ilnbuf)
         read(luint) buf, ibuf, nut
         do int = 1, nut
            call unpack2(ibuf(int),indi,indj)
            ovrlp(indi,indj) = buf(int)
            ovrlp(indj,indi) = buf(int)
         end do 
      end do

      close(luint,status = 'KEEP')

      return
      end

