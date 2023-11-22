      subroutine orbfile(nbas,coefs)
      implicit none

      integer ppp, muu
 
      integer nbas
 
      double precision coefs(nbas,nbas)

      open(800,file='orbsfile')
      do ppp = 1, nbas
         do muu = 1, nbas
            write(800,10) coefs(muu,ppp) 
         end do
      end do
      close(800)

 10   format(F20.15)

      return
      end
 
