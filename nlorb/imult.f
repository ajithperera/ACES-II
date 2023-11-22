      subroutine imult(trans,tdim1,tdim2,matin,idim1,idim2,
     & matout,odim1,odim2)
      implicit none

      integer tdim1, tdim2, trans(tdim1,tdim2), idim1, idim2,
     & matin(idim1,idim2), odim1, odim2, matout(odim1,odim2),
     & iii, jjj, kkk

      do iii = 1, odim1
         do jjj = 1, odim2
            matout(iii,jjj) = 0
            do kkk = 1, tdim2
               matout(iii,jjj) = matout(iii,jjj) +
     & trans(iii,kkk)*matin(kkk,jjj)
            end do
         end do
      end do
 
      return
      end

