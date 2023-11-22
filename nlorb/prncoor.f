      subroutine prncoor(natoms,crd,a2u,ach)
      implicit none

      integer natoms, iii

      integer ach(natoms)

      double precision crd(3,natoms), a2u

      write(*,10)

      do iii = 1, natoms
         write(*,99) iii, ach(iii), crd(1,iii)/a2u, crd(2,iii)/a2u,
     & crd(3,iii)/a2u
      end do

      write(*,*)
      write(*,*)

 99   format(2I6,3F12.7)
 10   format('Cartesian coordinates')

      return
      end
