      subroutine initize(coord,atchrg,mombf,centbf,natoms,nbas)
      implicit none

CCC hughes

C      integer iii

CCC hughes

      integer iintln, ifltln, iintfp, ialone, ibitwd

      integer natoms, nbas

      integer atchrg(natoms), mombf(nbas), centbf(nbas)

      double precision coord(3,natoms)

      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd

      call getrec(20,'JOBARC','COORD   ',natoms*3*iintfp,coord)
      call getrec(20,'JOBARC','CENTERBF',nbas,centbf)
      call getrec(20,'JOBARC','ATOMCHRG',natoms,atchrg)
      call getrec(20,'JOBARC','ANGMOMBF',nbas,mombf)

CCC hughes

C      do iii = 1, nbas
C         write(*,*) iii, centbf(iii)
C      end do

CCC hughes
      
      return
      end
