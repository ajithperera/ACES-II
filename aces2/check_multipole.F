      subroutine check_multipole(nao, nao2, scr1, scr2)
c     
      implicit none
#include "machsp.com"
c     
      integer maxprop, nao, nao2
      parameter (maxprop = 19)  ! 19 = 3 + 6 + 10
      double precision scr1(nao2), scr2(nao, nao), x
      integer nopert, ipert, i, j, icount
c     
      character*8 PERTSTR,STRING
      dimension pertstr(maxprop)
c     
      NOPERT=19
      PERTSTR(1)='DIPOLE_X'
      PERTSTR(2)='DIPOLE_Y'
      PERTSTR(3)='DIPOLE_Z'
c     
      PERTSTR(4)='QUAD_XX '
      PERTSTR(5)='QUAD_YY'
      PERTSTR(6)='QUAD_ZZ '
      PERTSTR(7)='QUAD_XY '
      PERTSTR(8)='QUAD_XZ '
      PERTSTR(9)='QUAD_YZ '
c     
      PERTSTR(10)='OCTUPXXX'
      PERTSTR(11)='OCTUPYYY'
      PERTSTR(12)='OCTUPZZZ'
      PERTSTR(13)='OCTUPXXY'
      PERTSTR(14)='OCTUPXXZ'
      PERTSTR(15)='OCTUPXYY'
      PERTSTR(16)='OCTUPYYZ'
      PERTSTR(17)='OCTUPXZZ'
      PERTSTR(18)='OCTUPYZZ'
      PERTSTR(19)='OCTUPXYZ'
c     
      do ipert = 1, nopert
         string = pertstr(ipert)
         call a2getrec(20,'JOBARC',string, nao2*iintfp, scr1)
         icount = 0
         do i = 1, nao
            do j = 1, i
               icount = icount + 1
               x = scr1(icount)
               scr2(i,j) = x
               scr2(j,i) = x
            enddo
         enddo
         write(6,*) ' Multipole ints ', string
         write(6,*)
         call output(scr2, 1, nao, 1, nao, nao, nao, 1)
      enddo
c     
      return
      end
