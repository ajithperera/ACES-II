c The entire operating procedure of symcor is explained in README.

      SUBROUTINE SYMCOR_MULDIM(ICORE,ICRSIZ)
      IMPLICIT INTEGER (A-Z)
      INTEGER ICORE(ICRSIZ), I0

C The dimension of IRRNM array is changed to 14 to accomodate the
C largest possible dimension. Note that in joda D_infh and C_infv
C are set to D_8h and C_8v. Ajith 08/2001

      DOUBLE PRECISION E, TRAD, PRD_ENRG_CHNG, TAU, TotEng,
     &                 reference_e, one
      CHARACTER*8 GRPSYM,IRRNM(14)
      character*8 labeleng(2),label(32),irrnmc(14)
      CHARACTER*4 DOIT
      LOGICAL FIRST_PASS,GRAD_EXIST

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
#include "control.com"
#include "parallel_aces.com"
#include <mrcc_heff.com>
#include "flags.h"
#include "jodaflags.com"


      DATA IZERO /0/

c   o initialize the free-space index in iCore
      I0=1
      ione = 1

C The values get assigned for PASS1 JOBARC record has very important
C function in joda when it is executing any calculation that involve
C numerical derivatives.
C
C When the execution transfers to here for the first time, joda
C has run one time as a single point calculation. There are no
C JOBARC records GRADIENT or PASS1. So GRAD_EXIST is false,
C PASS1 is 0 and GMTRYOPT is true or false depending on whether
C we are doing geometry optimizations or vibrational frequencies.
C
C Numerical geometry optimization - PASS1 is set to 1 until
C full set of gradients available. Then it is set to 0.
C
C Analytical gradient frequencies  - PASS1 is set to 1 in the
C very first time. It remains as 1 until the full Hessian is
C is available, then it is set to 0.
C
C Numerical gradient frequencies  - PASS1 get set to -1 in the
C first run. Before the next execution joda, this  get reset
C back to 1. When the full Hessian is available it get reset to 0.
C (-1 in the first call tell the code that it has to do energy
C points). 03/2003 Ajith Perera.

      CALL GETREC(-1,'JOBARC','PASS1',1,PASS1)
      CALL GETREC(0,'JOBARC','GRADIENT',LENGTH,ITMP)
      GRAD_EXIST = LENGTH.NE.-1
#ifdef _DEBUG_LVL0
      write(6,*) ' @symcor_muldim: pass1 upon entry ', pass1
#endif
c
c   o initialize /CONTROL/
      ENERONLY = IFLAGS2(h_IFLAGS2_grad_calc).EQ.2
      ROTPROJ  = IFLAGS (h_IFLAGS_fd_project).EQ.0
      GRADONLY = IFLAGS2(h_IFLAGS2_grad_calc).EQ.1

c   o detect dangerous point groups and avoid them like hell
      CALL GETCREC(20,'JOBARC','FULLPTGP',8,GRPSYM)
      IPOS = LINBLNK(GRPSYM(1:4))
      IF (PASS1.NE.1.AND.
     &    (GRPSYM(1:1).EQ.'C'.AND.GRPSYM(IPOS:IPOS).NE.'v'.OR.
     &     GRPSYM(1:3).EQ.'T h'.OR.
     &     GRPSYM(1:3).EQ.'T  '.OR.
     &     GRPSYM(1:1).EQ.'S'
     &    )
     &   ) THEN
         DOIT='COMP'
         WRITE(6,410)
 410     FORMAT(T3,'@SYMCOR: Full point group is dangerous! ',
     &             'Abelian subgroup will be used.')
         ITMP=1
         CALL PUTREC(20,'JOBARC','DANGERUS',1,ITMP)
      ELSE
         CALL GETREC(-1,'JOBARC','DANGERUS',1,ITMP)
         IF (IFLAGS(h_IFLAGS_fd_usegroup).EQ.0.AND.ITMP.EQ.0) THEN
            DOIT='FULL'
         ELSE
            DOIT='COMP'
            ITMP=0
         END IF
         IF (PASS1.NE.1) CALL PUTREC(20,'JOBARC','DANGERUS',1,ITMP)
      END IF

      CALL GETREC(20,'JOBARC','NATOMS',1,NATOM)

      IF (PASS1.EQ.-1) THEN
c
#ifdef _DEBUG_LVL0
      write(6,*) ' @symcor_muldim: entered pass1 = -1 block '
#endif
C
C In this block for purely numerical deriavtive calculations. Note that
C deriavtives of <A|DE/DQ|B> which is computed here only available 
C as finite differences. Ajith Perera, 04/2007. 
C
            CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,E)
C
C This call initilize the arrays necessary to compute aforementioned
C derivatives.
C
         CALL NUMFDS4_HEFF(ICORE, ICRSIZ)
C
         CALL PUTREC(20,'JOBARC','REFENERG',IINTFP,E)
         PASS1=1
         CALL PUTREC(20,'JOBARC','PASS1',1,PASS1)
      ELSE
         E = 0.d0
      END IF

      FIRST_PASS = ( (PASS1.EQ.0) .OR. (GRAD_EXIST.AND.GMTRYOPT) )

      IF (FIRST_PASS) THEN
#ifdef _DEBUG_LVL0
         write(6,*) '@symcor_muldim: Entered first_pass branch '
#endif

c      o create records that must survive truncating JOBARC
         CALL PUTREC(20,'JOBARC','REFENERG',IINTFP,E)
         CALL PUTREC(20,'JOBARC','LASTGEOM',1,IZERO)

c      - these three records are used by vee and vcceh
         CALL PUTREC(20,'JOBARC','PRINSPIN',1,IZERO)
         CALL PUTREC(20,'JOBARC','PRINFROM',IINTFP,E)
         CALL PUTREC(20,'JOBARC','PRININTO',IINTFP,E)
c
c      - these are used for "restarts"
         CALL PUTCREC(20,'JOBARC','PREVPTGP',4,'NONE')
         CALL PUTREC(0,'JOBARC','TGSSOCCA',8,0)
         CALL PUTREC(0,'JOBARC','TGSSOCCB',8,0)

c      o initialize the finite difference grid
         call init_fds4_heff(doit,natom,icore,icrsiz,nener,
     &                       IFLAGS(h_IFLAGS_print).GE.10)

         IF (.NOT.GMTRYOPT.AND.NENER.NE.0.AND.IRANK.EQ.0) THEN
c         o do the reference geom iff this is a vib freq calc with numerical
c           gradients and it is the root process
            PASS1=-1
         ELSE
            PASS1=1
         END IF
         CALL PUTREC(20,'JOBARC','PASS1',1,PASS1)
c
         write(6,*) ' GMTRYOPT ', GMTRYOPT
         write(6,*) ' NENER ', NENER
         write(6,*) ' irank ', irank
         write(6,*) ' Pass1 is set to ', pass1
c
cmn There are more things in heff.F that presumably should be written here ...
c
         call putrec(20,'JOBARC','PARENERG',
     $        iintfp,0.0d0)         
         call putrec(20,'JOBARC','REFER_E ',
     $        iintfp,reference_e)         
c     
c     create markers for dimensions of the matrices
c     
      ndim = 0
      call putrec(1, 'JOBARC', 'NDIABOCC', ione, ndim)
      call putrec(1, 'JOBARC', 'NDIABVRT', ione, ndim)
         do icase = 1, ncases
            ndim = 0
            call putrec(1, 'JOBARC', ndimheff(icase), ione, ndim)
            call putrec(1, 'JOBARC', namenvec(icase), ione, ndim)
         enddo
         call putrec(1, 'JOBARC','DIABFACT',15*iintfp, icore)
         call getrec(1, 'JOBARC', 'NATOMS', 1, natoms)
         call putrec(1, 'JOBARC','DIAB_GAP',natoms*iintfp, icore)
         call putrec(20, 'JOBARC', 'NUMBASIS', ione, ione)
c
cmn The record below is now taken care of by nextgeom ?  ... No something does not work
c
cmn         call putrec(1, 'JOBARC', 'JOBARCST', ione, ione)
cmn         call zeroja('JOBARCST')
c     
C
C The following record no matter how awkward it looks is essential. Note at
C this point we ran xjoda and there are certain records in JOBARC and then we
C add the "NEXTGEOM" record as a marker (initialized to junk at icore). Then we
C continue with vmol2ja and so on. In the next pass we are going  to fill the
C "NEXTGEOM" (see below the call to NEXTGEO) with proper data. The call to
C ACES_JA_TRUNCATE would delete all the records that comes after "NEXTGEOM"
C record. If we don't fixed the position of the "NEXTGEOM" record after the
C first JODA run, the records that are dependent on the symmetry get  written
C to the JOBARC before the "NEXTGEOM" appear and ACES_JA_TRUNCATE call would do
C nothing to erase them. Then When it times to run the lower symmtry points,
C the PUTREC is going to cry since we will be trying to change the record
C lengths!  Ajith Perera 08/2000.
C
         CALL PUTREC(20,'JOBARC','NEXTGEOM',3*NATOM*IINTFP,ICORE)
c
         write(6,*) ' @symcor_muldim: Process FDCALCTP, ipass ', pass1
         call getrec(20,'JOBARC','NUMPOINT',ione, npoint)
         CALL GETREC(20,'JOBARC','FDCALCTP',npoint,ICORE)
c
         call putrec(1, 'JOBARC', 'JOBARCST', ione, ione)
         call aces_ja_truncate('JOBARCST', 1)
         if (PASS1 .eq. -1) RETURN
      END IF
C
      ione = 1
      I000 = IONE
      call getrec(20,'JOBARC','NUMPOINT',ione, npoint)
      I010 = I000  + 3*NATOM*NPOINT*IINTFP
      I020 = I010 + NPOINT
      IF ((ICRSIZ-I020+I0) .LE. 0) CALL ACES_EXIT(1)

      write(6,*) ' Process FDCALCTP, ipass ', pass1
      CALL GETREC(20,'JOBARC','FDCALCTP',npoint,ICORE)
c
      write(6,*) ' Enter upd_fds4_heff '
      call upd_fds4_heff(NATOM, DOIT, IMORE, ICORE(I000), 
     &                  ICORE(I010), ICORE(I020), npoint,
     &                  (ICRSIZ-I020+I0)/IINTFP)
c
#ifdef _DEBUG_LVL0
      write(6,*) ' after @upd_fds4_heff, imore: ', imore
#endif
      IF (IMORE.EQ.1) THEN
c      o if there are more points, delete all the geometry-dependent records
c      o The record MARKZERO has been written lastly in upd_fds4_heff
         CALL ACES_JA_TRUNCATE('MARKZERO',1)
      ELSE
c
         if (iflags2(h_IFLAGS2_scan_normal) .eq. 1) then
            CALL GETREC(20,'JOBARC','NUMVIBRT',1, NMODE)
            iscan = 1
            iscr = iscan + nmode*iintfp
            write(6,*) ' Ready to process final results: scan_heff '
            call write_scan_heff(nmode, icore(iscan), icore(iscr), 
     &                  (ICRSIZ-Iscr+I0)/IINTFP)
c
         elseif (iflags2(h_IFLAGS2_geo_sample) .eq. 1) then
            iscr = 1
            write(6,*) ' Ready to process final results: geo_sample '
            call write_geo_heff(icore(iscr), 
     &                  (ICRSIZ-Iscr+I0)/IINTFP)
c
         elseif (iflags2(h_IFLAGS2_gen_quartic) .eq. 1) then
            iscr = 1
            CALL GETREC(20,'JOBARC','NUMVIBRT',1, NMODE)
            write(6,*) ' Ready to process final results: gen_quartic '
            call write_quartic(icore(iscr), nmode)
c
         elseif(iflags2(h_IFLAGS2_prep_vibron) .eq. 1) then
c
         write(6,*) ' Ready to process final results: prep_vibron '
c
         CALL GETREC(-1,'JOBARC',DOIT//'NIRX',1,NIRREPF)
         CALL GETREC(20,'JOBARC','NUMVIBRT',1, NMODE)

C
            Ifreq = IONE 
         inumirpf = ifreq + NMODE*iintfp
             iscr = inumirpf + nirrepf*iintfp
C
c         call write_eff(icore(ifreq), numvib, nirrepf,
c     $        icore(inumirpf), icore(iscr))
C
#ifdef _DEBUG_LVL0
         write(6,*) ' nirrepf, nmode, ifreq, iscr ',
     $            nirrepf, nmode, ifreq, iscr
#endif

         call write_pntheff(icore(iscr))
         write(6,*) ' Finished writing PNTHEFF file '
         call write_vibronic(icore(ifreq), nmode, icore(iscr))
         write(6,*) ' Finished writing Vibronic_Coupling file '
         call write_vibron_input(icore(ifreq), nmode, icore(iscr))
         write(6,*) ' Finished writing Vibron_Input file '

         else
            write(6,*) ' Unexpected case in symcor_muldim '
            call aces_exit(1)
         endif

         write(6,*) ' End in symcor_muldim '

         ione = 1
c      o signal the end of this findiff series
         CALL PUTREC(20,'JOBARC','LASTGEOM',1,ione)
         call putrec(1,'JOBARC','FNDFDONE',1,ione)
         call putrec(1,'JOBARC','JODADONE',1,ione)
         PASS1=0
         CALL PUTREC(20,'JOBARC','PASS1',1,PASS1)
      END IF

      RETURN
      END

