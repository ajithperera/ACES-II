      subroutine transform(trans,scr,matrix,nbastot)
      implicit none
 
      integer nbastot
      double precision matrix(nbastot,nbastot),
     & trans(nbastot,nbastot), scr(nbastot,nbastot)

      call dzeroarr(nbastot,nbastot,scr)
      call xgemm('N','N',nbastot,nbastot,nbastot,1.0D0,trans,
     & nbastot,matrix,nbastot,0.0D0,scr,nbastot)

      call dzeroarr(nbastot,nbastot,matrix)
      call xgemm('N','T',nbastot,nbastot,nbastot,1.0D0,scr,
     & nbastot,trans,nbastot,0.0D0,matrix,nbastot)
  
      return
      end
