      subroutine calc_tmom(dens, scr1, scr2, nao, nao2)
c     
      implicit none
c     
#include "machsp.com"
      
c     
      integer maxprop, nao, nao2
      parameter (maxprop = 19)  ! 19 = 3 + 6 + 10
      double precision scr1(nao2), scr2(nao, nao), x,
     $     xprop(maxprop), dens(nao*nao), sdot, y
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
c      PERTSTR(9)='2NDMO_YZ'
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
         call getrec(20,'JOBARC',string, nao2*iintfp, scr1)
         icount = 0
         do i = 1, nao
            do j = 1, i
               icount = icount + 1
               x = scr1(icount)
               scr2(i,j) = x
               scr2(j,i) = x
            enddo
         enddo
         x = SDOT(nao*nao, dens, 1, dens, 1)
         y = SDOT(nao*nao, scr2, 1, scr2, 1)
         xprop(ipert) = SDOT(nao*nao, dens, 1, scr2, 1)
c
         if (ipert .eq. 3 .and. .false.) then
            write(6,*) ' @calc_tmom, <Z> ', xprop(ipert)
            write(6,*) ' @calc_tmom, Trace Z-integrals ', y
c            call output(scr2, 1, nao, 1, nao, nao, nao, 1)
            write(6,*) ' @calc_tmom, Trace dens ', x
c            call output(dens, 1, nao, 1, nao, nao, nao, 1)
         endif
      enddo
c     
      call putrec(20,'JOBARC','MULTMOM2',
     $     maxprop*iintfp, xprop)
c     
      return
      end
