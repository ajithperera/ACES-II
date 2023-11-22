      subroutine setpts4_pes_scan(natom, nirrepf, type, label,
     &     ipttype, points,
     &     harm_freq, dscr, ndscr)
c     
c     this routine determines which points will be run and what type
c     of calculation (gradient or energy only) is most efficient for
c     each symmetry block
c     
C-----------------------------------------------------------------------
C     JDW  May/June 1996. Note on modifications of dimensioning.
C     
C     Memory allocation :
C     
C     It is assumed that there will never be more than (3*NATOM)**2
C     energy/gradient calculations. The memory allocation in the main
C     program for SETPTS and NXTGEO is based on this assumption.
C     Note that this bound is highly wasteful for GRADONLY calculations,
C     as the number of gradients will be much less. However, the code
C     currently reads record GRDPOINT even during ENERONLY jobs, hence
C     we need to make appropriate allowances.
C     
C     Dimensions allocated in main program (upper bounds).
C     
C     DSCR     : Double precision, 6*NATOM
C     POINTS  : Double precision, 9 * NATOM * NATOM * 3 * NATOM
C     (holds all sets of cartesian coordinates)
C     IPTTYPE : Integer, 9 * NATOM * NATOM
C     (a value for each point)
C-----------------------------------------------------------------------
C     
      IMPLICIT NONE

#include "flags.h"
#include "control.com"
#include "machsp.com"
#include "parallel_aces.com"
#include "jodaflags.com"
#include "mxatms.par"
#include "coord.com"

      INTEGER NATOM, NIRREPF, IPTTYPE, NDSCR

      CHARACTER*8 LABEL(NIRREPF)
      CHARACTER*4 TYPE, FULLREP, COMPREP, FNAMEIRP
      CHARACTER*1 a_ri
      CHARACTER*5 ANAME
      CHARACTER*80 FNAME
C     
      DOUBLE PRECISION POINTS, HARM_FREQ, DSCR
      INTEGER SKIPIR, NPTIRR, IJUNK, NMODES_PER_IRREP
      DIMENSION NMODES_PER_IRREP(14),HARM_FREQ(3*NATOM),
     &     POINTS(27*NATOM*NATOM*NATOM+9*NATOM*NATOM)
      DIMENSION IPTTYPE(9*NATOM*NATOM+3*NATOM), DSCR(NDSCR),
     &     nptirr(100),IJUNK(20), SKIPIR(20)
      LOGICAL PRINTQ, NORMAL_FDIF_PRSNT, match4, curvy_normal, print
      DOUBLE PRECISION TOL, ONE, ZILCH, FRED, REFERENCE_E,
     $     STPSIZ, XNORM, SCAN_STEP, STEP, FACT
      INTEGER ICOUNT, I, LFREE, LNCOR, LSCAN, LSCR, NDSCRLFT, lrefcoor,
     $     IUNIT, IOFFQ, IMODE, IRP, IOFFX, K, NVIB, NUMVIB, NMODE,
     $     MIN_SCAN, MAX_SCAN, IOFFP, NPOINT, IOFF, NSIZE, IONE, LCOOR,
     $     J, NTOT, iroot, ir, io, ne, ix, lsymq, lmass, lr,
     $     lq, ircurvy
      integer map_modes(200)
C     

C
C     
      DATA TOL   /1.D-8/
      DATA ONE   /1.0d0/
      DATA ZILCH /0.0/
C     
      fred = 0.0911355d0
      reference_e = 0.0d0

      PRINTQ=(IFLAGS(1).GE.10 .or. .true.)
      print = .true.
      NSIZE=3*NATOM
      ione = 1

      if (gmtryopt.or.gradonly.or.nprocs.eq.1) then
         iRoot = 0
      else
c     o the root process is doing the reference geometry in a vib freq calc
c     with numerical gradients
         iRoot = 1
      end if

c     o check for restrictions on FD irreps
      CALL GETREC(-1,'JOBARC','NFDIRREP',1,ICOUNT)
      IF (ICOUNT.NE.0) THEN
         CALL GETREC(-1,'JOBARC','FDIRREP ',ICOUNT,IJUNK)
         WRITE(6,501)
 501     FORMAT(T3,'@SETPTS: FCM evaluation limited ',
     &        'to the following symmetries :')
         WRITE(6,'((12I5))')(IJUNK(I),I=1,ICOUNT)
         DO I=1,20
            SKIPIR(I)=0
         END DO
         DO I=1,ICOUNT
            SKIPIR(IJUNK(I))=1
         END DO
         write(6,*) ' Array skipir '
         write(6,'((12i5))')(skipir(i),i=1,nirrepf)
      ELSE
         DO I=1,20
            SKIPIR(I)=1
         END DO
      END IF

      STPSIZ=DFLOAT(IFLAGS(57))*10.0D-5
      IF (PRINTQ) THEN
         WRITE(6,500)STPSIZ
 500     FORMAT(T3,'Step size will be ',F8.5,' amu**(1/2) * bohr.')
      END IF

      curvy_normal =  (iflags2(h_IFLAGS2_normal_curvy).eq.1) 
c       curvy_normal =  (iflags2(159).eq.1) 
c
      lFREE = 1
      
      lsymq = lfree
      lfree = lfree + nsize*nsize      
      lCOOR = lFREE
      lFREE = lFREE + NSIZE
      lrefCOOR = lFREE
      lFREE = lFREE + NSIZE
      lMASS = lFREE
      lFREE = lFREE + NATOM
      lNCOR = lFREE
      lFREE = lFREE + NSIZE*NSIZE 
      lSCAN = lFREE 
      lFREE = lFREE + NSIZE + 1
      lr = lFREE
      lFREE = LFREE + nsize
      lq = lFREE
      lFREE = LFREE + nsize
      ircurvy = lFREE
      lFREE = lFREE + nsize
      lSCR = lFREE
      lFREE = lFREE + NSIZE
      NDSCRLFT = NDSCR+1-lFREE
C     
      IF (NDSCRLFT.LT.0) THEN
         print *, '@SETPTS: Insufficient memory.'
         print *, '         need ',-ndscrlft*ifltln,' more bytes'
         call aces_exit(1)
      END IF
C     
C     The SCAN is built on top of normal modes of the parent state,
C     Read them from the NORMAL_FDIF file. 
c

c
      INQUIRE (FILE='NORMAL_FDIF', EXIST=NORMAL_FDIF_PRSNT)
C     
      If (NORMAL_FDIF_PRSNT) then
         write(6,*) '  Parent state normal coordinates are to be used'
      else
         write(6,*) ' At present Heff calculations require '
         write(6,*) ' NORMAL_FDIF file which contain '
         write(6,*) ' Normal modes of parent state in special format'
         write(6,*) ' This Run is Aborted.'
         call aces_exit(1)
      endif
c     
      CALL GETREC(20,'JOBARC','NUMVIBRT',1,numvib)
      CALL GETREC(20,'JOBARC',TYPE//'LABL',NIRREPF*IINTFP,LABEL)
c
      call rd_normal_fdif(natom, nirrepf,
     $     type, label, nmodes_per_irrep, skipir,
     $    harm_freq, lsymq, dscr, ndscr, map_modes,
     $     numvib)
      
c
      call putrec(20, 'JOBARC', 'FREQS_0 ', numvib*iintfp, harm_freq)
      call putrec(20, 'JOBARC', 'NMODES_F', nirrepf, nmodes_per_irrep)
c     
      CALL GETREC(20,'JOBARC','COORD   ',NSIZE*IINTFP,DSCR(lCOOR))
      CALL GETREC(20,'JOBARC','COORD   ',NSIZE*IINTFP,DSCR(lrefCOOR))
      CALL PUTREC(20,'JOBARC','REFCOORD',NSIZE*IINTFP,DSCR(lCOOR))

      CALL GETREC(20,'JOBARC','ATOMMASS',NATOM*IINTFP,DSCR(lMASS))
      CALL PUTREC(20,'JOBARC','REFMASS ',NATOM*IINTFP,DSCR(lMASS))

      CALL GETREC(20,'JOBARC','COORD   ',NSIZE*IINTFP,DSCR(lCOOR))

      IF (PRINTQ) THEN
         CALL GETREC(20,'JOBARC',TYPE//'LABL',NIRREPF*IINTFP,LABEL)
      END IF
      
      CALL GETREC(20,'JOBARC','NUMVIBRT',1,NMODE)

      if (curvy_normal) then
         write(6,*) ' overwrite normal modes with curvy_fdif '
         call rd_curvy_fdif(nsize, nsize-6, nmode, 
     $        harm_freq, dscr(lsymq), map_modes, dscr(lncor))
c
c we also need ZMAT coordinates in dscr(lcoor) ...
c
           Call Getrec(20, 'JOBARC', 'CORD_INT', NSIZE*IINTFP,
     $        dscr(lcoor))
           call sqush(dscr(lcoor), nsize)
      endif
      
      CALL ZERO(POINTS,27*NATOM*NATOM*NATOM+9*NATOM*NATOM)
C     
C     We don't get this far if SCAN_INFO file does not exsist. So we
C     can open it.

      iunit = 41
c     
c     SCAN_INFO should contain the following:
c     min_scan, max_scan, scan_step 
c     nmode
c     reference energy in a.u. 
c     displacement vector along reduced normal 
c     coordinates of parent state. One entry per line.
c
      write(6,*) ' Start reading SCAN_INFO file '
      open(iunit, file='SCAN_INFO')
      read(iunit, *) min_scan, max_scan, scan_step
      read(iunit,*) numvib
      if (numvib .ne. nmode) then
         write(6,*) ' # of modes on SCAN_INFO does not match NMODE'
         write(6,*) ' error in stepts4_scan_pes.F '
         call errex
      endif
      read(iunit,*) reference_e
      do i = 1, nmode
         read(iunit,*) dscr(lSCAN + i)
      enddo
c     
c     scr(i) contains the scan direction. It is supposed to be normalized.
c     
      xnorm = 0.0d0
      do i = 1, nmode
         xnorm = xnorm + dscr(i+lSCAN)**2
      enddo
C     
      if (abs(xnorm) .lt. 1.0d-2) then
         write(6,*) ' Input norm of scanning vector is too low'
         write(6,*) ' check SCAN_INFO file'
         call errex
      endif
      if (abs(xnorm - 1.0d0) .gt. 1.0d-10) then
         write(6,*) ' Input vector from SCAN_INFO will be normalized'
         xnorm = 1.0d0 / sqrt(xnorm)
         do i = 1, nmode
            dscr(i+lSCAN) = dscr(i+lSCAN) * xnorm
         enddo
      endif
c     
c     generate all displacements
c     
#ifdef _DEBUG_LVL0
      write(6,*) ' displacements from SCAN_INFO'
      call output(dscr(lSCAN+1), 1, nmode, 1, 1, nmode, 1, 1)
#endif
c     
      ioffp  = 1
      npoint = 0
      ioffq  = 1
C
      write(6,*) ' Harmonic frequencies in setpts4_pes_scan '
      call output(harm_freq, 1, 1, 1, nmode, 1, nmode, 1)
      write(6,*) ' Infamous Fred ', FRED
c
      do imode = 1, nmode
         step = dscr(imode+lSCAN) / (fred*sqrt(harm_freq(imode)))
#ifdef _DEBUG_LVL0
c         write(6,*) ' step along imode ', step, imode
c         write(6,*) ' normal mode coordinate'
c         call output(dscr(ioffq+(imode-1)*nsize), 1, 3,
c     $        1, natom, 3, natom, 1)
#endif
C     
         call SAXPY(nsize,
     $        step, dscr(ioffq+(imode-1)*nsize), 1,
     $        points(ioffp), 1)
      enddo
      if (print) then
      write(6,*) ' reference geometry'
      call output(dscr(lcoor), 1, 3, 1, natom, 3, natom ,1)
      write(6,*) ' basic unit displacement vector'
      call output(points(ioffp), 1, 3, 1, natom, 3, natom ,1)
      endif

c     
CSSS  call zero(points(ioffp), nsize)
      do i = min_scan, max_scan
         if (i .ne. 0) then
            call zero(points(ioffp), nsize)
c     
c     generate coordinates with + i * scan_step * Q 
c     
            ioffq = 1
            do imode = 1, nmode
               step = i * dscr(imode+lSCAN) * scan_step 
c     
c     normal modes need to be scaled to reduced normal modes
c     
               step = step / (fred*sqrt(harm_freq(imode)))
               call SAXPY(nsize,
     $              step, dscr(ioffq+(imode-1)*nsize), 1,
     $              points(ioffp), 1)
            enddo
            npoint = npoint + 1
c
#ifdef _DEBUG_LVL0
      write(6,*) ' write displacement vector for step', i
      call output(points(ioffp), 1, 3, 1, natom, 3, natom ,1)
#endif

            ioffp = ioffp + nsize
         endif
      enddo
c     
c make a copy of R and Q in COORD common blocks
c
      if (curvy_normal .or. .true.) then
      call dcopy(nsize, R, 1, dscr(lr), 1)
      call dcopy(nsize, Q, 1, dscr(lq), 1)
      endif
      ioff = 1
      do i = 1, npoint
         call SAXPY(nsize, one, dscr(lcoor), 1, points(ioff), 1)
         if (curvy_normal) then
            if (.true.) then
               write(6,*) ' Internal coordinates directly from S'
               call print_internal(points(ioff), nsize, nsize-6,
     $              natom, 6)
            endif
            call
     $        symcor_cart(dscr(lrefcoor),
     $           points(ioff), dscr(lscr), nsize, nsize-6)
         if (.true.) then
            write(6,*) 'internal coordinates for point ', i
            call dcopy(nsize, points(ioff), 1, Q, 1)
            call xtor(dscr(ircurvy),1)
         endif
         endif
         ioff = ioff + nsize
      enddo
      ntot = npoint
c
c restore R,Q to initial geomtry; probably not necessary ?  No! it is.
c
      if (curvy_normal .or. .true.) then
      call dcopy(nsize, dscr(lr), 1, R, 1)
      call dcopy(nsize, dscr(lq), 1, Q, 1)
      endif
c     
      write(6,1000)npoint
 1000 format(t3,' total number of single-point energy calculations:',i5)
      if(npoint.eq.0)then
         write(6,1005)
 1005    format(t3,'@setpts-f, there are no vibrational modes!')
         call aces_exit(1)
      endif

      write(6,*) ' Parallel points genererated ',
     $     iRank, nprocs, iroot
      iOff = 0
      do iR = 0, nProcs-1
         call paces_batch_stat(iR,nProcs,iRoot,nPoint,iO,nE)
         write(6,*) ' iO, nE ', iO, nE
         do iX = 1+iO, nE+iO
            iPtType(iOff+iX) = 1+iR
         end do
      end do
      iRoot = mod(iRoot+nPoint,nProcs)

      if (print) then
      do i = 1, npoint
         write(6,*) ' ipttype ', ipttype(i), i
      enddo
      write(6,*) ' Set of displacements'
      ioff = 1
      do i = 1, npoint
         write(6,*) ' Geometry ', i
         call output(points(ioff), 1, 3, 1, natom, 3, natom, 1)
         ioff = ioff + nsize
      enddo
      endif
c     
      CALL PUTREC(20,'JOBARC','FDCOORDS',NTOT*NSIZE*IINTFP,POINTS)
      CALL PUTREC(20,'JOBARC','FDCALCTP',NTOT,IPTTYPE)
      CALL PUTREC(20,'JOBARC','NUMPOINT',IONE,NPOINT)

      CALL IZERO(NPTIRR,NIRREPF)
      CALL PUTREC(20,'JOBARC','NPTIRREP',NIRREPF,NPTIRR)
c     
      return
      end
