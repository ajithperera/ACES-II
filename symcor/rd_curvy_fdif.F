      subroutine rd_curvy_fdif(NRX1,NX6M, nmode, Sfreq, Snormal,
     $      map_modes, SCR)
c     
c     read in symmetry adapted internal coordinates from ACES, using ZMAT definition
c
c     nrx1 = 3 * natoms
c     nx6m = 3 * natoms - 6
c     natoms includes ghost atoms
c
      implicit none

      integer NX6M, NRX0, NX6M0
      integer NRX, NRX6M, nmode, NRX1
      double precision Sfreq(NRX1), Snormal(NX6M,NRX1),
     $     SCR(NRX1,NRX1)
      character*10 VARNAM
      integer map_modes(200)

      integer iunit, itop, ibot, j, icoord, i, irot
      logical yesno
      character*80 string
      logical print
c     
      print = .false.
c     
      iunit = 70
      inquire(file='CURVY_FDIF', exist = yesno)
      if (.not. yesno) then
         write(6,*) ' Expect to find file CURVY_FDIF '
         call errex
      endif
      call zero(Snormal, nx6m*nrx1)
      call zero(Sfreq, nrx1)
c     
      OPEN(UNIT=70,FILE='CURVY_FDIF')
      REWIND(iunit)

      IBOT=1
      read(iunit,*) string
 3117 FORMAT(T20,'   Normal modes in internal coordinates ')
      read(iunit,*) NRX0, NX6m0
      if (NX6m0 .ne. NX6M) then
         write(6,*) ' inconsistency dimensions curvy_fdif, ZMAT '
         write(6,*) ' @RD_CURVY_FDIF : nx6m, NX6M0',
     $        NX6M, NX6M0
         call errex
      endif
      NRX = NRX0
c
 3114 ITOP=MIN(IBOT+4,NRX)
      READ(iunit,*) string
 3116 FORMAT(75('-'))
      READ(iunit,*)(Sfreq(j),J=IBOT,ITOP)
      DO 3115 ICOORD=1,NX6M
         READ(iunit,*)VARNAM,
     &        (Snormal(ICOORD,J),J=IBOT,ITOP)
 3115 CONTINUE
      IBOT=ITOP+1
      IF(IBOT.LT.NRX)GOTO 3114
      READ(iunit,*)
      close(70, status='KEEP')
c
      if (print) then
          write(6,*) ' Frequencies of Internal normal modes '
         call output(Sfreq, 1, 1, 1, nrx, 1, nrx, 1)
         write(6,*) ' Internal coordinates directly read from file '
         call output(Snormal, 1, nx6m, 1, nrx, nrx1, nrx, 1)
      endif
c     
c     resort coordinates to represent true normal modes 
c
      irot = 0
 100  irot = irot + 1

      if (abs(Sfreq(irot)) .lt. 1.0d-3
     $     .and. irot .lt. nrx) goto 100
      irot = irot - 1
c
      do i = irot + 1, nrx
         j = i - irot
         Sfreq(j) = Sfreq(i)
         call DCOPY(nx6m, Snormal(1,i), 1, Snormal(1, j), 1)
      enddo
c
      if (nmode .ne. nrx - irot) then
         write(6,*) ' Something unexpected rd_curvy_fdif '
         write(6,*) ' nmode, nrx, irot, (nrx-irot) ',
     $        nmode, nrx, irot, nrx-irot
         write(6,*) ' Additional low frequencies are discarded ?'
         write(6,*) ' Better check everything is fine '
c         call errex
      endif
c
c copy normal modes in proper order (same as in normal_fdif)
c
      call dzero(SCR, NRX1*NRX1)
      do i = 1, nmode
         call DCOPY(NX6M, Snormal(1,i), 1, SCR(1,i), 1)
      enddo
      do i = 1, nmode
         j = map_modes(i)
         call DCOPY(nx6m, SCR(1,j), 1, Snormal(1,i), 1)
      enddo
      call DCOPY(nmode, Sfreq, 1, SCR, 1)
      do i = 1, nmode
         j = map_modes(i)
         Sfreq(i) = SCR(j,1)
      enddo
c      
      if (print) then
         write(6,*) ' Resorted frequencies and normal modes '
         write(6,*) ' Frequencies of Internal normal modes '
         call output(Sfreq, 1, 1, 1, nmode, 1, nmode, 1)
         write(6,*) ' Internal coordinates processed from file '
         call output(Snormal, 1, nx6m, 1, nmode, nrx1, nmode, 1)
      endif
c      
c
c print info on frequencies
c
      write(6,*)
      write(6,*) ' Vibrations read from curvy_fdif '
      write(6,*) ' index new   index old    frequency  '
      write(6,*) ' ----------------------------------- ' 
      do i = 1, nmode
         write(6,880) i, map_modes(i), Sfreq(i)
      enddo
 880  format(i6, 4x, i6, 4x, F14.6)
c
c in the calling subroutine Snormal has leading dimension NRX1, as SCR here. Copy!
c
      call dzero(SCR, NRX1*NRX1)
      do i = 1, nmode
         call DCOPY(NX6M, Snormal(1,i), 1, SCR(1,i), 1)
      enddo
c
         call DCOPY(NRX1*nmode, SCR, 1, Snormal, 1)
c
c also print out resorted normal modes, using proper leading dimension
c
      write(6,*) ' Resorted Curvilinear normal modes '
      call output(Snormal, 1, nx6m, 1, nmode, nrx1, nmode, 1)
c
      return
      end
