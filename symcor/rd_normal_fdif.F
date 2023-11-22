      subroutine rd_normal_fdif(natom, nirrepf,
     $     type, label, nmodes_per_irrep, skipir,
     $    harm_freq, lsymq, dscr, ndscr, map_modes,
     $     numvib)
c
      implicit none
      integer natom, nirrepf, ndscr, lsymq
      CHARACTER*8 LABEL(NIRREPF)
      CHARACTER*4 TYPE, FULLREP, COMPREP, FNAMEIRP
      CHARACTER*1 a_ri
      CHARACTER*5 ANAME
      CHARACTER*80 FNAME
      INTEGER NMODES_PER_IRREP(14)
      Double Precision fact, dscr(ndscr),HARM_FREQ(3*NATOM)
C     
      integer iunit, imode, irp, i, skipir(20), lscr, ioffx, 
     $     map_modes(200), nvib, ioffq, numvib, ioff, k, j, nsize
      logical match4
c
      iunit = 70
      open(unit=iunit,file='NORMAL_FDIF',status='OLD')
      call zero(harm_freq, 3*natom)
      ioffq = lsymq
      nsize = 3*natom
      lscr = lsymq + nsize*numvib
c
      imode = 0
      do 50 irp = 1, nirrepf
         rewind(iunit)
         nmodes_per_irrep(irp) = 0
         if (skipir(irp) .eq. 0) goto 50
C     
         do i = 1, 4
            fnameirp(i:i) = label(irp)(i:i)
         enddo
c     

         write(6,*) ' Processing modes of irrep', irp, fnameirp
         do 100 i = 1, numvib
c            write(6,*) ' mode i, total ', i, numvib
c            call flush(6)
            read(iunit,2111) fullrep, comprep, fact, a_ri
            write(6,2111) fullrep, comprep, fact, a_ri
            ioffx = lSCR 
            do k = 1, natom
c               write(6,*) ' atom k ', k
c               call flush(6)
               read(iunit, 2115) aname, (dscr(ioffx+j), j=0,2)
c               write(6,*) ' done read ', k
c               call flush(6)
               ioffx = ioffx + 3
            enddo
C     
 2111 FORMAT(T3,A4,T12,A4,T20,F14.2,T40,A1)
 2115 FORMAT(T3,A5,T10,F14.8,T30,F14.8,T50,F14.8) 
c     
            if (match4(fullrep,fnameirp,4)
     $     .and. abs(fact) .gt. 30.0D0
     $           .and. a_ri .ne. 'i') then
               nmodes_per_irrep(irp) = nmodes_per_irrep(irp) + 1
               imode = imode + 1
               harm_freq(imode) = fact
               map_modes(imode) = i
               call SCOPY(nsize, dscr(lSCR), 1, dscr(ioffq), 1)
               ioffq = ioffq + nsize
            endif
c     
 100     continue
 50   continue
c
      close(unit=iunit,status='KEEP')
      write(6,*) ' Modes included per irrep'
      nvib = 0
      do irp = 1, nirrepf
         nvib = nvib + nmodes_per_irrep(irp)
         write(6,*) irp, nmodes_per_irrep(irp)
      enddo
      if (nvib .ne. numvib) then
         write(6,*) ' Number of included vibrations has changed'
         write(6,*) ' New number of vibrations ', nvib
         numvib = nvib
         CALL PUTREC(20,'JOBARC','NUMVIBRT',1,numvib)
      endif
c
c write map_modes to jobarc
c
      CALL PUTREC(20,'JOBARC','MAP_MODE',numvib,map_modes)
c
c print info on frequencies
c
      write(6,*)
      write(6,*) ' Vibrations read from normal_fdif '
      write(6,*) ' index new   index old    frequency  '
      write(6,*) ' ----------------------------------- ' 
      do i = 1, nvib
         write(6,880) i, map_modes(i), harm_freq(i)
      enddo
 880  format(i6, 4x, i6, 4x, F14.6)
c
c also print out resorted normal modes
c
      write(6,*) ' Resorted Normal modes '
      call output(dscr(lsymq), 1, nsize, 1, nvib, nsize, nvib, 1)
c
      return
      end
