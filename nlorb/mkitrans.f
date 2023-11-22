      subroutine mkitrans(trans,nbas,nbastot,scr,row,numrow)
      implicit none

      integer nbas, nbastot, trans(nbastot,nbas),
     & scr(nbas,nbas), numrow, row(numrow)

      call izeroarr(scr,nbas)
      call icpident(scr,nbas)
      
      call delrow(row,numrow,scr,trans,nbas)

      return
      end
