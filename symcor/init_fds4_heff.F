
c     This routine performs all the symmetry analysis needed to generate the grid
c     of displacements. It is essentially a method that uses and updates JOBARC
c     records. NENER is returned solely to tell symcor that there are points to
c     start working on (or not).

c     INPUT
c     char*4  DOIT
c     integer NATOM
c     integer ICRSIZ
c     logical BPRINT

c     OUTPUT
c     integer ICORE(ICRSIZ)
c     integer NENER

c     RECORDS
c     get DOIT//'ORDR'
c     get DOIT//'NORB'
c     get DOIT//'NIRX'
c     get 'COMPORDR'
c     get 'COMPNORB'
c     get 'COMPNIRX'



      subroutine init_fds4_heff(doit,natom,icore,icrsiz,nener,bprint)
      implicit none

#include "control.com"
#include "parallel_aces.com"
#include "flags.h"
#include "jodaflags.com"

      character*4 doit
      integer natom, icore(*), icrsiz, nener, nmodes_per_irrep(14),
     &     length, npoints_per_irrep(14) 
      logical bprint, scan_info_prsnt, scan_pes, geo_sample, gen_quartic
      

      integer iorder, iorderf, iorderc
      integer norbit, norbitf, norbitc
      integer nirrep, nirrepf, nirrepc
      integer I000, I010, I020, I030, I040, I050, I060, I070, I080, I090
      integer I100, I110, I120, I130, I140, I150, I160, I170, I180, I190
      integer I200, I0, i
      character*8 irrnm(32), label(32), irrnmc(32)
      character*80 fname

#include "control.com" /* for ROTPROJ */
#include "machsp.com"

      CALL GETREC(20,'JOBARC',DOIT//'ORDR',1,IORDERF)
      CALL GETREC(20,'JOBARC',DOIT//'NORB',1,NORBITF)
      CALL GETREC(20,'JOBARC','COMPORDR',1,IORDERC)
      CALL GETREC(20,'JOBARC','COMPNORB',1,NORBITC)
      IORDER=MAX(IORDERF,IORDERC)
      NORBIT=MAX(NORBITF,NORBITC)
cmn
      nener = 1
      I0 = 1

      I000=I0
      I010=I000 + IINTFP*IORDER*9
      I020=I010 + IINTFP*IORDER*IORDER*9
      I030=I020 + IINTFP*IORDER*IORDER
      I040=I030 + IINTFP*IORDER*IORDER
      I050=I040 + 3*IORDER
      I060=I050 + IORDER*IORDER
      I080=I060 + IORDER + IAND(IORDER,1)
      I090=I080 + 6*IORDER
      I100=I090 + IINTFP*3*IORDER

      CALL CHRTABLE(IORDERF,ICORE(I000),ICORE(I010),ICORE(I020),
     &     ICORE(I030),ICORE(I040),ICORE(I050),ICORE(I060),
     &     IRRNM,ICORE(I080),DOIT,ICORE(I090))
      CALL CHRTABLE(IORDERC,ICORE(I000),ICORE(I010),ICORE(I020),
     &     ICORE(I030),ICORE(I040),ICORE(I050),ICORE(I060),
     &     IRRNMC,ICORE(I080),'COMP',ICORE(I090))

c     
c     write the names of the irreps to a file
c     
      open(unit=70, file='IRPNAMES')
      rewind(70)
      write(70,*) ' Names of irreps full symmetry group'
      write(70,*) iorderf
      do i = 1, iorderf
         write(70,*) irrnm(i)
      enddo
      write(70, *)
      write(70,*) ' Names of irreps Computational symmetry group'
      write(70,*) iorderc
      do i = 1, iorderc
         write(70,*) irrnmc(i)
      enddo
      close(70)
c     

      CALL GETREC(-1,'JOBARC',DOIT//'NIRX',1,NIRREPF)
      CALL GETREC(-1,'JOBARC','COMPNIRX',1,NIRREPC)

      call putrec(20, 'JOBARC', 'NAMEIRPF', nirrepf, irrnm)
      call putrec(20, 'JOBARC', 'NAMEIRPC', nirrepc, irrnmc)
c     
      NIRREP=MAX(NIRREPF,NIRREPC)

      I010=I000 + IINTFP*NATOM*NATOM*9
      I020=I010 + IINTFP*NATOM*NATOM*9
      I030=I020 + IINTFP*NATOM*NATOM*9
      I040=I030 + IINTFP*MAX(3*NATOM,NIRREPF*IORDER)
      I050=I040 + IINTFP*MAX(6*NATOM,9*IORDER)
      I060=I050 + IINTFP*MAX(6*NATOM,NIRREPF)
      I060=I050
      I070=I060 + IINTFP*9*NATOM*NATOM*3*NATOM
      I080=I070 + MAX(MAX(2,IORDER)*NATOM,9*NATOM*NATOM)
      I090=I080 + NATOM
      I100=I090 + NATOM + IAND(NATOM,1)
      I110=I100 + IINTFP*NATOM*NATOM*MAX(18,IORDER)
      I120=I110 + IINTFP*NATOM*NATOM*18
      I130=I120 + IINTFP*NATOM*NATOM*18

      CALL VIBINF(NATOM,NIRREPF,IORDERF,ICORE(I000),
     &     ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &     LABEL,ICORE(I060),ICORE(I070),ICORE(I080),
     &     ICORE(I090),ICORE(I100),DOIT)
      CALL SYMADQ(NATOM,NIRREPC,IORDERC,NORBITC,ICORE(I000),
     &     ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &     LABEL,ICORE(I060),ICORE(I070),ICORE(I080),
     &     ICORE(I090),ICORE(I100),ICORE(I110),ICORE(I120),
     &     ICORE(I130),'COMP',1)

      CALL TRAPRJ(NATOM,ICORE(I000),ICORE(I010),ICORE(I020),
     &     ICORE(I030))

      IF (ROTPROJ) THEN
         CALL ROTPRJ(NATOM,ICORE(I000),ICORE(I010),ICORE(I020),
     &        ICORE(I030))
      END IF

      CALL SCHMIDT(NATOM,ICORE(I010),ICORE(I020))

      CALL COLLECT(NATOM,NIRREPC,ICORE(I000),
     &     ICORE(I010),ICORE(I020),ICORE(I030),LABEL,
     &     'COMP')

      IF (BPRINT) THEN
         CALL PRTCOORD(NATOM,NIRREPC,0,ICORE(I000),ICORE(I010),
     &        LABEL,'COMP',.FALSE.)
      END IF

      CALL SYMADQ(NATOM,NIRREPF,IORDERF,NORBITF,ICORE(I000),
     &     ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &     LABEL,ICORE(I060),ICORE(I070),ICORE(I080),
     &     ICORE(I090),ICORE(I100),ICORE(I110),ICORE(I120),
     &     ICORE(I130),DOIT,2)
      CALL PRDEGEN(NATOM,NIRREPF,IORDERF,NORBITF,ICORE(I000),
     &     ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &     LABEL,ICORE(I060),ICORE(I070),ICORE(I080),
     &     ICORE(I090),ICORE(I100),DOIT)

      IF (BPRINT) THEN
         CALL PRTCOORD(NATOM,NIRREPF,NIRREPC,ICORE(I000),ICORE(I010),
     &        LABEL,DOIT,.TRUE.)
      END IF

      I010=I000
      I020=I010 + 9*NATOM*NATOM + 3*NATOM
      I030=I020 + IINTFP*(9*NATOM*NATOM*3*NATOM + 9*NATOM*NATOM)
      I200=I030 + 3*NATOM * IINTFP

      scan_pes = iflags2(h_IFLAGS2_scan_normal) .eq. 1
      geo_sample = iflags2(h_IFLAGS2_geo_sample) .eq. 1
      gen_quartic = iflags2(h_IFLAGS2_gen_quartic) .eq. 1
      if (scan_pes) then
c     
c     We will create the points for a PES scan. Read the displacement vector
c     info from the file "SCAN_INFO"
c     
         INQUIRE (FILE='SCAN_INFO', EXIST=SCAN_INFO_PRSNT)
         if (SCAN_INFO_PRSNT) then
            write(6,*) '  Read in Scan info'
         else
            write(6,*) ' At present Scan calculations require '
            write(6,*) ' Scan_info file which contains '
            write(6,*) ' the info on displacements '
            write(6,*) ' See symcor/setpts4_pes_scan on format'
            write(6,*) ' Please provide this file and try again '
            call aces_exit(1)
         endif

         CALL SETPTS4_PES_SCAN(NATOM,NIRREPF,DOIT,IRRNM,
     &        ICORE(I010),ICORE(I020),ICORE(I030),
     &        ICORE(I200),
     &        (ICRSIZ-I200+I0)/IINTFP)
C     
      elseif (geo_sample) then
c     
c     We will create the points for a PES sampling. Read the displacement vector
c     info from the file "GEO_SCAN"
c     
         INQUIRE (FILE='GEO_SAMPLE', EXIST=SCAN_INFO_PRSNT)
         if (SCAN_INFO_PRSNT) then
            write(6,*) '  Read in GEO_SAMPLE'
         else
            write(6,*) ' At present Scan calculations require '
            write(6,*) ' GEO_SAMPLE file which contains '
            write(6,*) ' the info on displacements '
            write(6,*) ' See symcor/setpts4_pes_scan on format'
            write(6,*) ' Please provide this file and try again '
            call aces_exit(1)
         endif

         CALL SETPTS4_PES_SAMPLE(NATOM,NIRREPF,DOIT,IRRNM,
     &        ICORE(I010),ICORE(I020),ICORE(I030),
     &        ICORE(I200),
     &        (ICRSIZ-I200+I0)/IINTFP)
C     
      elseif (gen_quartic) then
c     
c     We will create the points for a Quartic force field. Read INFO
c     info from the file "GEN_QUARTIC"
c     
         INQUIRE (FILE='GEN_QUARTIC', EXIST=SCAN_INFO_PRSNT)
         if (SCAN_INFO_PRSNT) then
            write(6,*) '  Read in GEN_QUARTIC'
         else
            write(6,*) ' At present gen_quartic calculations require '
            write(6,*) ' GEN_QUARTIC file which contains '
            write(6,*) ' the info on displacements '
            write(6,*) ' See symcor/setpts4_quartic on format file'
            write(6,*) ' Please provide this file and try again '
            call aces_exit(1)
         endif

         CALL SETPTS4_QUARTIC(NATOM,NIRREPF,DOIT,IRRNM,
     &        ICORE(I010),ICORE(I020),ICORE(I030),
     &        ICORE(I200),
     &        (ICRSIZ-I200+I0)/IINTFP)
C     
      else
         CALL SETPTS4_HEFF(NATOM,NIRREPF,DOIT,IRRNM,
     &        ICORE(I010),ICORE(I020),ICORE(I030),
     &        ICORE(I200),
     &        (ICRSIZ-I200+I0)/IINTFP)
      endif

      return
c     end subroutine init_fds4_heff
      end

