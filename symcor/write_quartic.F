      subroutine write_quartic(freq, numvib)
c     
      implicit none
#include <machsp.com>
c     
      integer ione, npoint, nmode, iunit, nquartic, numvib,
     $     isample, jmode0, imode0, i, map_modes(200), i0
      DOUBLE PRECISION Q1(5), Q2(2,15), delta,
     $     scr(30), freq(numvib)
c     
      DATA Q1 / -2.0, -1.0, 0.0, 1.0, 2.0 /
      DATA Q2 /  0.0,  0.0,
     $     -2.0,  0.0,
     $     -1.0,  0.0,
     $     1.0,  0.0,
     $     2.0,  0.0,
     $     0.0, -2.0,
     $     0.0, -1.0,
     $     0.0,  1.0,
     $     0.0,  2.0,
     $    -2.0,  2.0,
     $     2.0,  1.0,
     $     1.0,  2.0,
     $     -2.0, -2.0,
     $     -2.0, -1.0,
     $     -1.0, -2.0 /
c     
      ione = 1
      CALL GETREC(20,'JOBARC','NUMPOINT',IONE,npoint)
      CALL GETREC(20,'JOBARC','NUMVIBRT',1,NMODE)
      call getrec(1, 'JOBARC', 'FREQS_0 ', numvib*iintfp,
     $     freq)
      CALL GETREC(20,'JOBARC','MAP_MODE',nmode,map_modes)
C     
C     We don't get this far if GEN_QUARTIC file does not exist. So we
C     can open it.

      iunit = 41
c     
c     GEN_quartic should contain the following:
c     nquartic
c     delta
c     nmode (the total number of normal modes (for checkin purposes)
c     imode, jmode (jmode if nquartic = 2)
c     
      write(6,*) ' Start reading GEN_QUARTIC file '
      open(iunit, file='GEN_QUARTIC')
      rewind(iunit)
      read(iunit,*)  delta
      read(iunit, *) nquartic
      write(6,*) ' nquartic ', nquartic
c     
      if (.true.) then
         read(iunit,*) numvib
         if (numvib .ne. nmode) then
            write(6,*) ' # of modes on GEN_Quartic does not match NMODE'
            write(6,*) ' error in write_quartic.F '
            call errex
         endif
      endif
      if (nquartic .eq. 1) then
         read(iunit,*) imode0
         isample = 5
         jmode0 = 0
      else
         read(iunit,*) imode0, jmode0
         isample = 15
      endif
c     
      close(iunit,status='KEEP')
c     
      if (imode0 .ne. 0) then
         do i = 1, nmode
            if (map_modes(i) .eq. imode0) i0 = i
         enddo
         imode0 = i0
      endif
      if (jmode0 .ne. 0) then
         do i = 1, nmode
            if (map_modes(i) .eq. jmode0) i0 = i
         enddo
         jmode0 = i0
      endif
c
      if (isample .ne. npoint) then
         write(6,*) ' Something wrong in write_quartic :',
     $        npoint, isample
         call errex
      endif
c     
c     read in the calculated energies:
c     
      call getrec(20,'JOBARC','PNTENERG',
     $     npoint*iintfp, scr)
c     
      iunit = 42
      open(42, file='E_QUARTIC')
      write(iunit,*) nquartic
      if (nquartic .eq. 1) then
         write(iunit, 885) imode0, freq(imode0) / 8065.54092884153d0,
     $        freq(imode0)
         elseif (nquartic .eq. 2) then
         write(iunit, 886) imode0, freq(imode0) / 8065.54092884153d0,
     $           freq(imode0), jmode0,
     $           freq(jmode0) / 8065.54092884153d0, freq(jmode0)
      endif
 885  format(i6, F14.8)
 886  format(i6, F14.8, F14.4, 6x, i6, F14.8, F14.4)
      write(iunit,*)
      write(iunit,*) isample      
      if (nquartic .eq. 1) then
         do i = 1, isample
            write(iunit, 880) i, Q1(i)*delta,
     $           (scr(i) -scr(1))*27.2113956
         enddo
      else
         do i = 1, isample
            write(iunit, 881) i, Q2(1,i)*delta, Q2(2,i)*delta,
     $           (scr(i) -scr(1))*27.2113956
         enddo
      endif
      close(iunit, status='KEEP')
c     
      if (nquartic .eq. 1) then
         write(6,80) imode0
 80      format(' Quartic displacements along normal mode ', i6)
      else
         write(6,81) imode0, jmode0
 81      format(' Quartic displacements along normal modes ', 2i6)
      endif
c     
      write(6,*)
      write(6,*) ' Total energies for set of quartic displacements '
      write(6,*) ' ----------------------------------------------- '
      write(6,*)
      if (nquartic .eq. 1) then
         write(6, 82)
      else
         write(6,83)
      endif
 82   format(' point # ', 4x, ' Q1 ', 16x, ' Relative energy (eV) ')
 83   format(' point # ', 4x, ' Q1 ', 16x, ' Q2 ', 16x,
     $     ' Relative energy (eV) ')
      if (nquartic .eq. 1) then
         do i = 1, isample
            write(6, 880) i, Q1(i)*delta,
     $           (scr(i) -scr(1))*27.2113956
         enddo
      else
         do i = 1, isample
            write(6, 881) i, Q2(1,i)*delta, Q2(2,i)*delta,
     $           (scr(i) -scr(1))*27.2113956
         enddo
      endif
      write(6,*)
      write(6,*) ' ----------------------------------------------- '
      write(6,*)
      
c     
 880  format(i6, F20.12, F25.12)
 881  format(i6, 2F20.12, F25.12)
c     
      return
      end
      
