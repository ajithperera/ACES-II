      subroutine int_occupy(evalc,int_occ)

c Determines the occupancy needed for 
c the the integrator if it is a KS SCF calculation.

#include <aces.h>
      implicit none
#include <iuhf.com>
#include <symm2_ks.com>
#include <mol.com>

      integer int_occ(nbastot,2)
      M_REAL
     &    evalc(nbastot,2)
      
      integer ispin,i,l,k,m

      call callstack_push('INT_OCCUPY')
     
C prkash
        call nocc2occup

      do 100 ispin=1,1+iuhf
        l=0
        call izero(int_occ(1,ispin),nbastot)
         do 90 i=1,nirrep
          if (i.gt.1) then
            l=l+numbasir(i-1)
          end if
          k=occup(i,ispin)
          if (k.eq.0) go to 90
            do m=1,k
              int_occ(l+m,ispin)=1
            end do

  90    continue

        call piksr2(nbastot,evalc(1,ispin),int_occ(1,ispin))
  100  continue

      call callstack_pop
      return
      end
