      subroutine setpts4_heff(natom, nirrepf, type, label,
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
C     SCR     : Double precision, NDSCR
C     POINTS  : Double precision, (9 * NATOM * NATOM + 3*natom) * 3 * NATOM
C     (holds all sets of cartesian coordinates)
C     IPTTYPE : Integer, 9 * NATOM * NATOM + 3*natom
C     (a value for each point)
C-----------------------------------------------------------------------
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*8 LABEL(NIRREPF)
      CHARACTER*4 TYPE, FULLREP, COMPREP, FNAMEIRP
      CHARACTER*1 A_RI
      CHARACTER*5 ANAME
      CHARACTER*80 FNAME
C     
      DIMENSION POINTS(27*NATOM*NATOM*NATOM
     &     +9*NATOM*NATOM)
      DIMENSION IPTTYPE(9*NATOM*NATOM+3*NATOM), DSCR(NDSCR),
     &     NPOINTS_PER_IRREP(14), HARm_FREQ(3*NATOM), 
     &     NMODES_PER_IRREP(14)
      INTEGER IJUNK(20), ICOUNT, SKIPIR(20), IDEGEN(100), 
     &     NPTIRR(100), ioffirp(14)
      LOGICAL PRINTQ, NORMAL_FDIF_PRSNT, DEGENERATE, print,
     $     match4, curvy_normal, reduced
      integer map_modes(200)
C     
#include "flags.h"
#include "control.com"
#include "machsp.com"
#include "parallel_aces.com"
#include "jodaflags.com"
#include "mxatms.par"
#include "coord.com"

C     
      DATA TOL   /1.D-8/
      DATA ONE   /1.0d0/
      DATA ZILCH /0.0/
C     
      fred = 0.0911355d0
      reference_e = 0.0d0
      print = .true.
      curvy_normal =  (iflags2(h_IFLAGS2_normal_curvy).eq.1) 

      PRINTQ=(IFLAGS(1).GE.10)
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
         WRITE(6,'(12I5)')(IJUNK(I),I=1,ICOUNT)
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
      
c     
c     monitor step size. if negative we use displacments in dimensionless reduced coordinates.
      if (IFLAGS(57) .lt. 0) then
         write(6,*)  ' Use displacements along ',
     $        'dimensionless reduced coordinates.'
         write(6,*) ' Use Infamous Fred for scaling ', fred
         reduced = .true.
         stpsiz = - DFLOAT(IFLAGS(57)) / 1000.0d0
         write(6,511) stpsiz
 511     format(T3,' Elementary stepsize ', F8.3, ' Dimensionless ')
      else
         STPSIZ=DFLOAT(IFLAGS(57))*10.0D-5
         write(6,*) ' Use Normal mode displacements '
         WRITE(6,500)STPSIZ
 500     FORMAT(T3,'Step size will be ',F8.5,' amu**(1/2) * bohr.')
      endif

      lFREE = 1
      lsymq = lfree
      lfree = lfree + nsize*nsize
      lCOOR = lFREE
      lFREE = lFREE + NSIZE
      lREFCOOR = lFREE
      lFREE = lFREE + NSIZE
      lMASS = lFREE 
      lFREE = lFREE + NATOM
      lr = lFREE
      lFREE = LFREE + nsize
      lq = lFREE
      lFREE = LFREE + nsize
      lSCR  = lFREE
      lFREE = lFREE + 3*NATOM
      NDSCRLFT = NDSCR+1-lFREE
C     
      IF (NDSCRLFT.LT.0) THEN
         print *, '@SETPTS4_HEFF: Insufficient memory.'
         print *, '         need ',-ndscrlft*ifltln,' more bytes'
         call aces_exit(1)
      END IF
C     
C     The grid is built on top of normal modes of the parent state.
C     Read them from the NORMAL_FDIF file.
C     
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
C     
      CALL GETREC(20,'JOBARC','NUMVIBRT',1,numvib)
      CALL GETREC(20,'JOBARC',TYPE//'LABL',NIRREPF*IINTFP,LABEL)
      
c     
      call rd_normal_fdif(natom, nirrepf,
     $     type, label, nmodes_per_irrep, skipir,
     $     harm_freq, lsymq, dscr, ndscr, map_modes,
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
      
      CALL GETREC(20,'JOBARC',TYPE//'DEGN',NIRREPF,IDEGEN)
      
      CALL ZERO(POINTS,27*NATOM*NATOM*NATOM+9*NATOM*NATOM)
      CALL IZERO(NPTIRR,NIRREPF)

c     o set up vector of reciprocal square roots of atomic masses
      DO IOFF=0,NATOM-1
         X=SQRT(DSCR(lMASS+IOFF))
         IF (X.LT.TOL) THEN
            DSCR(lMASS+IOFF)=0.d0
         ELSE
            DSCR(lMASS+IOFF)=1.d0/X
         END IF
      END DO
c     
      write(6,*)
      write(6,*) ' initial coordinates '
      call output(dscr(lcoor), 1, 3, 1, natom, 3, natom, 1)
      write(6,*)
c     
      if (curvy_normal) then
         write(6,*) ' overwrite normal modes with curvy_fdif '
         call rd_curvy_fdif(nsize, nsize-6, numvib, 
     $        harm_freq, dscr(lsymq), map_modes, dscr(lfree))
         write(6,*) ' Done reading curvy_fdif '
         call flush(6)
c     
c     we also need ZMAT coordinates in dscr(lcoor) ...
c     
         Call Getrec(20, 'JOBARC', 'CORD_INT', NSIZE*IINTFP,
     $        dscr(lcoor))
         call sqush(dscr(lcoor), nsize)
      endif
c     
      call zero(points,27*natom*natom*natom+9*natom*natom)
      CALL GETREC(20,'JOBARC','NUMVIBRT',1, numvib) 
c     
      ioffirp(1) = 0
      do irp = 2, nirrepf
         ioffirp(irp) = ioffirp(irp-1) + nmodes_per_irrep(irp-1)
      enddo

      ioffp = 1
      npoint = 0
c     
c     generate coordinates with +Q 
c     
      ioffq = 1
      do irp = 1, nirrepf
         do i=1, nmodes_per_irrep(irp)
            iq = ioffirp(irp) + i
            if (print) then
               write(6,*) ' Plus displacement '
               write(6,*) ' Normal mode, irp, i ', irp, i
               call output(dscr(ioffq+(i-1)*nsize),
     $              1, 1, 1, nsize, 1, nsize, 1)
            endif
            call zero(points(ioffp), nsize)
            if (reduced) then
               step = stpsiz / (fred*sqrt(harm_freq(iq)))
            else
               step = stpsiz
            endif
            write(6,882), iq, harm_freq(iq), step
 882        format(' normal mode ', i6, ' frequency ', F12.4,
     $           ' step ', E12.4)
            call SAXPY(nsize,
     $           step, dscr(ioffq+(i-1)*nsize), 1,
     $           points(ioffp), 1)
            ioffp=ioffp+nsize
            npoint = npoint + 1
            iq = ioffirp(irp) + i
            jq = 0
            write(6,880) npoint, iq, jq
         enddo
         ioffq = ioffq + nmodes_per_irrep(irp)*nsize
         npoints_per_irrep(irp) = npoint
      enddo
c     
 880  format(' point ', i4, ' modea ', i4, ' modeb ', i4)
c     
c     generate coordinates with -Q
c     
c     if we have a degenerate point group we run over all + / - displacements
c     if the point group is non-degenerate then we only do - displacements 
c     for the symmetric representation
c     
c     first determine if this is a degenerate subgroup
c     
#ifdef _NOSKIP
C     This block of code generate deg_map it is not used (as far as I can tell).
C     SO SKIP. If the info is needed the use IDGEN array which contain the same
C     info and directly read from JOBARC. 05/2007, Ajith Perera.
C     
      degenerate = .false.
c     
      nirr = 0
      do imode = 1, numvib
         deg_map(imode) = 0
      enddo
      do imode = 1 ,numvib
         if (deg_map(imode) .eq. 0) then
            nirr = nirr + 1
            deg_map(imode) = nirr
            f1 = freq(imode)
            do jmode = imode + 1, numvib
               if (abs(f1 - freq(jmode))
     $              .lt. 1.0d-5) then
                  deg_map(jmode) = nirr
                  write(6,*) ' degenerate mode found', imode, jmode
                  write(6,*) ' frequencies ', f1, freq(jmode)
               endif
            enddo
         endif
      enddo
c     
      degenerate = nirr .ne. numvib
c     
      if (degenerate) then
         nirrepf0 = nirrepf
      else
         nirrepf0 = nirrepf
      endif
c     
      write(6,*) ' @setpts_eff : degenerate pointgroup ? ',
     $     degenerate
      write(6,*) ' @setpts_eff : determine nirrepf0 ', nirrepf0
c     
#endif
C     
      Degenerate = .FALSE.
      DO Irrep = 1, Nirrepf
         Degenerate = (Idegen(Irrep) .GT. 1) 
      Enddo
c     
      write(6,*) ' npoint after +q displacements ', npoint
C     
      NIrrepf0 = NIrrepf
C     
      if (iflags2(h_IFLAGS2_grid_vibron) .ne. 1) then
c     
         ioffq = 1
         do irp = 1, nirrepf0
            do i=1, nmodes_per_irrep(irp)
               iq = ioffirp(irp) + i
               if (print) then
                  write(6,*) ' Minus displacement '
                  write(6,*) ' Normal mode, irp, i ', irp, i
                  call output(dscr(ioffq+(i-1)*nsize),
     $                 1, 1, 1, nsize, 1, nsize, 1)
               endif
               call zero(points(ioffp), nsize)
               if (reduced) then
                  step = stpsiz / (fred*sqrt(harm_freq(iq)))
               else
                  step = stpsiz
               endif
               call SAXPY(nsize,
     $              -step, dscr(ioffq+(i-1)*nsize), 1,
     $              points(ioffp), 1)
               ioffp=ioffp+nsize
               npoint = npoint + 1
               iq = -(ioffirp(irp) + i)
               jq = 0
               write(6,880) npoint, iq, jq
            enddo
            ioffq = ioffq + nmodes_per_irrep(irp)*nsize
            npoints_per_irrep(irp) =  npoints_per_irrep(irp) + npoint
         enddo
      endif
      
      write(6,*) ' npoint after -q displacements ', npoint
c     
      if (iflags2(h_IFLAGS2_grid_vibron) .eq. 8 .or. 
     $     iflags2(h_IFLAGS2_grid_vibron) .eq. 9) then
c     
c     generate coordinates +Qi, +Qj with same symmetry for degenerate irreps (off-diagonal only)
c     
         ioffq = 1
         do irp = 1, nirrepf
            do i = 1, nmodes_per_irrep(irp)
               do j = i+1, nmodes_per_irrep(irp)
                  iq = ioffirp(irp) + i
                  jq = ioffirp(irp) + j
                  if (abs(harm_freq(iq) - harm_freq(jq))
     $                 .lt. 1.0d-5) then
                     call zero(points(ioffp), nsize)
                     if (reduced) then
                        step = stpsiz / (fred*sqrt(harm_freq(iq)))
                     else
                        step = stpsiz
                     endif
                     call SAXPY(nsize,
     $                    step, dscr(ioffq+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     call SAXPY(nsize,
     $                    step, dscr(ioffq+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     write(6,880) npoint, iq, jq
                  endif
               enddo
            enddo
            ioffq = ioffq + nsize*nmodes_per_irrep(irp)
            npoints_per_irrep(irp) = npoints_per_irrep(irp) + npoint 
         enddo
c     
      endif
c     
      if (iflags2(h_IFLAGS2_grid_vibron) .eq. 3 .or. 
     $     iflags2(h_IFLAGS2_grid_vibron) .eq. 4 .or.
     $     iflags2(h_IFLAGS2_grid_vibron) .eq. 6 .or. 
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 7 .or.
     $     iflags2(h_IFLAGS2_grid_vibron) .eq. 10 .or. 
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 11) then
c     
c     generate coordinates +Qi, +Qj with same symmetry (off-diagonals only)
c     
         ioffq = 1
         do irp = 1, nirrepf
            do i = 1, nmodes_per_irrep(irp)
               do j = i+1, nmodes_per_irrep(irp)
                  iq = ioffirp(irp) + i
                  jq = ioffirp(irp) + j
                  call zero(points(ioffp), nsize)
                  if (reduced) then
                     step = stpsiz / (fred*sqrt(harm_freq(jq)))
                  else
                     step = stpsiz
                  endif
                  call SAXPY(nsize,
     $                 step, dscr(ioffq+(j-1)*nsize), 1,
     $                 points(ioffp), 1)
                  if (reduced) then
                     step = stpsiz / (fred*sqrt(harm_freq(iq)))
                  else
                     step = stpsiz
                  endif
                  call SAXPY(nsize,
     $                 step, dscr(ioffq+(i-1)*nsize), 1,
     $                 points(ioffp), 1)
                  ioffp = ioffp + nsize
                  npoint = npoint + 1
                  iq = ioffirp(irp) + i
                  jq = ioffirp(irp) + j
                  write(6,880) npoint, iq, jq
               enddo
            enddo
            ioffq = ioffq + nsize*nmodes_per_irrep(irp)
            npoints_per_irrep(irp) = npoints_per_irrep(irp) + npoint 
         enddo
c     
      endif
c     
      write(6,*) ' npoint after symmetric +qi +qj displacements ',
     $     npoint
c     
      if (iflags2(h_IFLAGS2_grid_vibron) .eq. 4 .or. 
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 7 .or.
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 11) then
c     
c     generate coordinates +Qi, +Qj with different symmetry (off-diagonals only)
c     
         ioffi = 1
         do irpi = 1, nirrepf
            ioffj = 1
            do irpj = 1, irpi - 1
               do i = 1, nmodes_per_irrep(irpi)
                  do j = 1, nmodes_per_irrep(irpj)
                     iq = ioffirp(irpi) + i
                     jq = ioffirp(irpj) + j
                     call zero(points(ioffp), nsize)
                     if (reduced) then
                        step = stpsiz / (fred*sqrt(harm_freq(jq)))
                     else
                        step = stpsiz
                     endif
                     call SAXPY(nsize,
     $                    step, dscr(ioffj+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     if (reduced) then
                        step = stpsiz / (fred*sqrt(harm_freq(iq)))
                     else
                        step = stpsiz
                     endif
                     call SAXPY(nsize,
     $                    step, dscr(ioffi+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     iq = ioffirp(irpi) + i
                     jq = ioffirp(irpj) + j
                     write(6,880) npoint, iq, jq
                  enddo
               enddo
               ioffj = ioffj + nsize*nmodes_per_irrep(irpj)
            enddo
            ioffi = ioffi + nsize*nmodes_per_irrep(irpi)
            npoints_per_irrep(irp) = npoints_per_irrep(irp) + npoint 
         enddo
      endif
c     
      write(6,*) ' npoint after asymmetric +qi +qj displacements ',
     $     npoint

      if (iflags2(h_IFLAGS2_grid_vibron) .eq. 8 .or. 
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 9) then
c     
c     generate coordinates -Qi, -Qj with same symmetry for degenerate irreps (off-diagonals only)
c     
         ioffq = 1
         do irp = 1, nirrepf
            do i = 1, nmodes_per_irrep(irp)
               do j = i+1, nmodes_per_irrep(irp)
                  iq = (ioffirp(irp) + i)
                  jq = (ioffirp(irp) + j)
                  if (abs(harm_freq(iq) - harm_freq(jq))
     $                 .lt. 1.0d-5) then
                     call zero(points(ioffp), nsize)
                     if (reduced) then
                        step = stpsiz / (fred*sqrt(harm_freq(iq)))
                     else
                        step = stpsiz
                     endif
                     call SAXPY(nsize,
     $                    -step, dscr(ioffq+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     call SAXPY(nsize,
     $                    -step, dscr(ioffq+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     write(6,880) npoint, -iq, -jq
                  endif
               enddo
            enddo
            ioffq = ioffq + nsize*nmodes_per_irrep(irp)
            npoints_per_irrep(irp) = npoints_per_irrep(irp) + npoint
         enddo
c     
      endif                     ! points .eq. 8 or 9


      if (iflags2(h_IFLAGS2_grid_vibron) .eq. 3 .or. 
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 4 .or.
     $     iflags2(h_IFLAGS2_grid_vibron) .eq. 6 .or. 
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 7 .or.
     $     iflags2(h_IFLAGS2_grid_vibron) .eq. 10 .or. 
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 11) then
c     
c     generate coordinates -Qi, -Qj with same symmetry (off-diagonals only)
c     
         ioffq = 1
         do irp = 1, nirrepf
            do i = 1, nmodes_per_irrep(irp)
               do j = i+1, nmodes_per_irrep(irp)
                  iq = ioffirp(irp) + i
                  jq = ioffirp(irp) + j
                  call zero(points(ioffp), nsize)
                  if (reduced) then
                     step = stpsiz / (fred*sqrt(harm_freq(jq)))
                  else
                     step = stpsiz
                  endif
                  call SAXPY(nsize,
     $                 -step, dscr(ioffq+(j-1)*nsize), 1,
     $                 points(ioffp), 1)
                  if (reduced) then
                     step = stpsiz / (fred*sqrt(harm_freq(iq)))
                  else
                     step = stpsiz
                  endif
                  call SAXPY(nsize,
     $                 -step, dscr(ioffq+(i-1)*nsize), 1,
     $                 points(ioffp), 1)
                  ioffp = ioffp + nsize
                  npoint = npoint + 1
                  iq = -(ioffirp(irp) + i)
                  jq = -(ioffirp(irp) + j)
                  write(6,880) npoint, iq, jq
               enddo
            enddo
            ioffq = ioffq + nsize*nmodes_per_irrep(irp)
            npoints_per_irrep(irp) = npoints_per_irrep(irp) + npoint
         enddo
c     
      endif                     ! points .eq. 3
c     
      write(6,*) ' npoint after symmetric -qi -qj displacements ',
     $     npoint
C     
      if (iflags2(h_IFLAGS2_grid_vibron) .eq. 4 .or. 
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 7 .or.
     &     iflags2(h_IFLAGS2_grid_vibron) .eq. 11) then
c     
c     generate coordinates -Qi, -Qj with different symmetry (off-diagonals only)
c     
         ioffi = 1
         do irpi = 1, nirrepf
            ioffj = 1
            do irpj = 1, irpi - 1
               do i = 1, nmodes_per_irrep(irpi)
                  do j = 1, nmodes_per_irrep(irpj)
                     iq = ioffirp(irpi) + i
                     jq = ioffirp(irpj) + j
                     call zero(points(ioffp), nsize)
                     if (reduced) then
                        step = stpsiz / (fred*sqrt(harm_freq(jq)))
                     else
                        step = stpsiz
                     endif
                     call SAXPY(nsize,
     $                    -step, dscr(ioffj+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     if (reduced) then
                        step = stpsiz / (fred*sqrt(harm_freq(iq)))
                     else
                        step = stpsiz
                     endif
                     call SAXPY(nsize,
     $                    -step, dscr(ioffi+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     iq = -(ioffirp(irpi) + i)
                     jq = -(ioffirp(irpj) + j)
                     write(6,880) npoint, iq, jq
                  enddo
               enddo
               ioffj = ioffj + nsize*nmodes_per_irrep(irpj)
            enddo
            ioffi = ioffi + nsize*nmodes_per_irrep(irpi)
            npoints_per_irrep(irp) = npoints_per_irrep(irp) + npoint
         enddo
c     
      endif                     ! points .eq. 4
c     
      write(6,*) ' npoint after asymmetric -qi -qj displacements ',
     $     npoint
c     
c     
      if ((iflags2(h_IFLAGS2_grid_vibron) .ge. 5)
     $     .and. (iflags2(h_IFLAGS2_grid_vibron) .ne. 8)) then
c     
c     generate coordinates with + 2 Q   (for cubic and quartic couplings)
c     
         ioffq = 1
         do irp = 1, nirrepf
            do i=1, nmodes_per_irrep(irp)
               iq = ioffirp(irp) + i
               call zero(points(ioffp), nsize)
               if (reduced) then
                  step = stpsiz / (fred*sqrt(harm_freq(iq)))
               else
                  step = stpsiz
               endif
               call SAXPY(nsize,
     $              2.0d0*step, dscr(ioffq+(i-1)*nsize), 1,
     $              points(ioffp), 1)
               ioffp=ioffp+nsize
               npoint = npoint + 1
               iq = ioffirp(irp) + i
               jq = ioffirp(irp) + i
               write(6,880) npoint, iq, jq
            enddo
            ioffq = ioffq + nmodes_per_irrep(irp)*nsize
            npoints_per_irrep(irp) = npoints_per_irrep(irp) + npoint
         enddo
c     
         write(6,*) ' npoint after +2 qi displacements ', npoint
c     generate coordinates with - 2 Q   (for cubic and quartic couplings)
c     
         if (degenerate) then
            nirrepf0 = nirrepf
         else
            nirrepf0 = nirrepf
         endif
c     
         ioffq = 1
         do irp = 1, nirrepf0
            do i=1, nmodes_per_irrep(irp)
               iq = ioffirp(irp) + i
               call zero(points(ioffp), nsize)
               if (reduced) then
                  step = stpsiz / (fred*sqrt(harm_freq(iq)))
               else
                  step = stpsiz
               endif
               call SAXPY(nsize,
     $              -2.0d0*step, dscr(ioffq+(i-1)*nsize), 1,
     $              points(ioffp), 1)
               ioffp=ioffp+nsize
               npoint = npoint + 1
               iq = -(ioffirp(irp) + i)
               jq = -(ioffirp(irp) + i)
               write(6,880) npoint, iq, jq
            enddo
            ioffq = ioffq + nmodes_per_irrep(irp)*nsize
            npoints_per_irrep(irp) = npoints_per_irrep(irp) + npoint
         enddo
c     
         write(6,*) ' npoint after -2 qi displacements ', npoint
      endif
c     
      if (iflags2(h_IFLAGS2_grid_vibron) .ge. 9) then
c     
c     generate coordinates aQi, bQj with same symmetry for degenerate irreps 
c     Used to calculate off-quartic couplings for degenerate irreps.
c     
         ioffq = 1
         do irp = 1, nirrepf
            do i = 1, nmodes_per_irrep(irp)
               do j = i+1, nmodes_per_irrep(irp)
                  iq = (ioffirp(irp) + i)
                  jq = (ioffirp(irp) + j)
                  if (abs(harm_freq(iq) - harm_freq(jq))
     $                 .lt. 1.0d-5) then
                     if (reduced) then
                        step = stpsiz / (fred*sqrt(harm_freq(iq)))
                     else
                        step = stpsiz
                     endif
c     -2 +2
                     call zero(points(ioffp), nsize)
                     call SAXPY(nsize,
     $                    -2*step, dscr(ioffq+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     call SAXPY(nsize,
     $                    2*step, dscr(ioffq+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     write(6,881) npoint, -2, iq, 2, jq
c     2 1
                     call zero(points(ioffp), nsize)
                     call SAXPY(nsize,
     $                    2*step, dscr(ioffq+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     call SAXPY(nsize,
     $                    1*step, dscr(ioffq+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     write(6,881) npoint, 2, iq, 1, jq
c     1 2
                     call zero(points(ioffp), nsize)
                     call SAXPY(nsize,
     $                    1*step, dscr(ioffq+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     call SAXPY(nsize,
     $                    2*step, dscr(ioffq+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     write(6,881) npoint, 1, iq, 2, jq
c     -2 -2
                     call zero(points(ioffp), nsize)
                     call SAXPY(nsize,
     $                    -2*step, dscr(ioffq+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     call SAXPY(nsize,
     $                    -2*step, dscr(ioffq+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     write(6,881) npoint, -2, iq, -2, jq
c     -2 -1
                     call zero(points(ioffp), nsize)
                     call SAXPY(nsize,
     $                    -2*step, dscr(ioffq+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     call SAXPY(nsize,
     $                    -1*step, dscr(ioffq+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     write(6,881) npoint, -2,iq, -1,jq
c     
c     -1 -2
                     call zero(points(ioffp), nsize)
                     call SAXPY(nsize,
     $                    -1*step, dscr(ioffq+(i-1)*nsize), 1,
     $                    points(ioffp), 1)
                     call SAXPY(nsize,
     $                    -2*step, dscr(ioffq+(j-1)*nsize), 1,
     $                    points(ioffp), 1)
                     ioffp = ioffp + nsize
                     npoint = npoint + 1
                     write(6,881) npoint, -1,iq, -2,jq
c     
                  endif
               enddo
            enddo
 881        format(' npoint = ', i4, ' stepa ', i3, ' modea ', i4,
     $           ' stepb ', i3, ' modeb ', i4)
     $           
            ioffq = ioffq + nsize*nmodes_per_irrep(irp)
            npoints_per_irrep(irp) = npoints_per_irrep(irp) + npoint
         enddo
c     
      endif                     ! points .eq. 9

c     
      if (npoint .gt. 9*natom*natom + 3*natom) then
         write(6,*) ' Too many points generated in setpts_eff'
         write(6,*) ' maximum ', 9*natom*natom + 3*natom
         write(6,*) ' generated points ', npoint
         call aces_exit(1)
      endif

C     
C     WRITE GEOMETRIES TO BE USED IN CALCULATION AND CALCULATION TYPES
C     TO THE JOBARC FILE
C     
      if (curvy_normal .or. .true.) then
         call dcopy(nsize, R, 1, dscr(lr), 1)
         call dcopy(nsize, Q, 1, dscr(lq), 1)
      endif
c     
      ioff = 1
      do i = 1, npoint
         call SAXPY(nsize, one, dscr(lcoor), 1, points(ioff), 1)
         if (curvy_normal) then
            call symcor_cart(dscr(lrefcoor),
     $           points(ioff), dscr(lscr), nsize, nsize-6)
         else
            call mn_align_geom(natoms, atmass,
     $           points(ioff), dscr(lrefcoor), dscr(lscr),
     $           ndscrlft)
         endif
         ioff = ioff + nsize
      enddo
      ntot = npoint
c     
c     restore R,Q to initial geomtry; probably not necessary ?
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
c     
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
c     
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

      write(6,*) ' @setpts4_heff.F, ntot, npoint ', ntot, npoint
c     
      CALL PUTREC(20,'JOBARC','FDCOORDS',NTOT*NSIZE*IINTFP,POINTS)
      CALL PUTREC(20,'JOBARC','FDCALCTP',NTOT,IPTTYPE)
      CALL PUTREC(20,'JOBARC','NUMPOINT',ione,NPOINT)
c     
      write(6,*) ' Check getrec fdcalctp'
      CALL GETREC(20,'JOBARC','FDCALCTP',NTOT,IPTTYPE)

      CALL IZERO(NPTIRR,NIRREPF)
      CALL PUTREC(20,'JOBARC','NPTIRREP',NIRREPF,NPTIRR)
c     
      return
      end
