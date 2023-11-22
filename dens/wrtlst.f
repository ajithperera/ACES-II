      subroutine wrtlst(Z,length)
      double precision Z(2)
c
c     write (6,10) length 
c10   format (/5x,'====  matrix in vdens with len=',i7)
      write (6,100) (Z(I),I=1,length)
 100  format (2x,10f7.4)
c
      return
      end
