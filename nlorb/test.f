      subroutine test(trans,transinv,nbastot)
      implicit none
 
      integer nbastot, iii, jjj, junk

      double precision trans(nbastot,nbastot),
     & transinv(nbastot,nbastot)

      junk = 0

      write(*,*) 
      write(*,*) 
 
      do iii = 1, nbastot
         do jjj = 1, nbastot
            if (trans(jjj,iii).ne.transinv(iii,jjj)) then
               write(*,10)
               stop
            end if
            if (trans(jjj,iii).eq.transinv(iii,jjj)) then
               junk = 1
            end if
         end do
      end do
      
      write(*,*) 
      write(*,*) 

 10   format('Reordering Transformation is Not Orthogonal')
 
      return
      end

 
