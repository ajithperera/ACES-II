      subroutine fillcoefvir(loc,can,nocc,nvir,nbastot)
      implicit none

      integer ppp, bbb

      integer nocc, nvir, nbastot

      double precision loc(nbastot,nbastot), can(nbastot,nbastot) 

      do ppp = 1, nbastot
         do bbb = 1, nvir
            loc(ppp,bbb+nocc) = 0.0D0
         end do
      end do

      do ppp = 1, nbastot
         do bbb = 1, nvir
            loc(ppp,bbb+nocc) = can(ppp,bbb+nocc) 
         end do
      end do

      return
      end

