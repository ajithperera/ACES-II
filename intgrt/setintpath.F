
          subroutine setintpath(relsiz,bslrd)

c  this routine determines the relative size of the atoms
c  partition using bragg-slater radii
c  equal partition
c  default is bragg-slater way
c  this routine is diferant from old routine as it does not have dynamic
C  partitioning as it involves knowledge of C and hence use of values from
c  primvec


       implicit none
#include <aces.h>
#include <intgrtflags.com>
#include <grid.com>
#include <mol.com>

      INTEGER iatm,jatm
      M_REAL
     &  relsiz(natoms,natoms),bslrd(natoms)

      call callstack_push('SETINTPATH')
      if(int_partpoly.eq.int_partpoly_bsrad)then

        do iatm=1,natoms
          do jatm=1,natoms
            if (iatm.ne.jatm) then
              relsiz(iatm,jatm)=bslrd(iatm)/(bslrd(iatm)+bslrd(jatm))
            end if
          end do
        end do

c Set relative size of atoms to 0.5
      else
        do iatm=1,natoms
          do jatm=1,natoms
            if (iatm.ne.jatm) then
               relsiz(iatm,jatm)=0.5d0
            end if
          end do
        end do

      end if

      call callstack_pop
      return
      end

