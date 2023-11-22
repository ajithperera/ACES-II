      subroutine transform2(trans,scr,matrix,nbastot)
      implicit none
 
      integer nbastot
      double precision matrix(nbastot,nbastot),
     & trans(nbastot,nbastot), scr(nbastot,nbastot)

      call dzeroarr(nbastot,nbastot,scr)
      call xgemm('N','N',nbastot,nbastot,nbastot,1.0D0,trans,
     & nbastot,matrix,nbastot,0.0D0,scr,nbastot)

      call dzeroarr(nbastot,nbastot,matrix)
      call dcopyarr(scr,matrix,nbastot)

      return
      end
