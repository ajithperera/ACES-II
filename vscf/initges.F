
      SUBROUTINE INITGES(ONEH,   XFORM,
     &                   DENS,   SCR1,
     &                   SCR2,   EVAL,
     &                   SCRA,   SCRB,
     &                   SCRTMP, SCRN,
     &                   CSO,    CZAO,
     &                   ZAOSO,  CNEW,
     &                   IMAGE,  ANGTYP,
     &                   OCCNUM, REFCOORD,
     &                   COORD,  ATMASS,
     &                   LDIM1,  LDIM2,
     &                   NBAS,NBASX,NATOMS,IUHF,IOLDOCC,
     &                   DCORE,  MAXDCOR)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C A complete overhaul by Ajith Perera to improve the readability
C and the flexiblity. 04/2005.
C-----------------------------------------------------------------------
      INTEGER          LDIM1,LDIM2,NBAS,NBASX,NATOMS,IUHF,MAXDCOR
      INTEGER          IMAGE(NBASX),ANGTYP(NBASX),OCCNUM(NBAS*(IUHF+1)),
     &                 IOLDOCC(16)
      DOUBLE PRECISION ONEH(LDIM1),XFORM(LDIM1),DENS((IUHF+1)*LDIM1),
     &                 SCR1(NBAS*NBAS),SCR2(LDIM2),EVAL((IUHF+1)*NBAS),
     &                 SCRA(LDIM2),SCRB(LDIM2),SCRTMP(LDIM2),
     &                 SCRN(NBAS*NBAS),CSO(NBAS*NBAS*(IUHF+1)),
     &                 CZAO(NBASX*NBAS,(IUHF+1)),ZAOSO(NBAS*NBASX),
     &                 CNEW(NBASX),REFCOORD(3,NATOMS),COORD(3,NATOMS),
     &                 ATMASS(NATOMS),DCORE(MAXDCOR)
C
      DOUBLE PRECISION ZILCH,ZTEST,DNRM2,DTMP,DDOT,Eng
      INTEGER          IONE,IKICK,NOCC,SWAP,LOCK,IPRTGS,ISTOP1,ISTOP2,
     &                 READMO,WRITMO,IUHFRHF,LUGSS,I,J,N,ICNT,
     &                 IFLAG45,IUNIT,PERLINE,NOCC_IRREP(8),INDX2,
     &                 Length,IScr,IOS
C
      LOGICAL ROHF,AOBASMOS_FILE,OLDMOS_FILE,Geom_Opt,Findf_vib,
     &        BRUKMOS_EXIST,ROHFMO,GSSOPT,GSSALW,hfstab,
     &        GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR,Occ_In_Zmt,
     &        Brueckner,G_FLE_READ,OLDAOMOS_FILE,AUTOMATE,
     &        AODENS_FILE, QRHF
C
      CHARACTER*11 GUESS_TYPE
C-----------------------------------------------------------------------
#include "machsp.com"
#include "molcas.com"
#include "symm2.com"
C
#include "flags.h"
#include "flags.com"
#include "flags2.com"
      COMMON /POPUL/  NOCC(16)
      COMMON /GSCOMA/ GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      COMMON /GSCOMB/ SWAP(4,8,2),LOCK(8,2),IPRTGS(8,2),ISTOP1,ISTOP2,
     &                READMO,WRITMO,IUHFRHF,LUGSS
      COMMON /GSCOMC/ ROHFMO
C-----------------------------------------------------------------------
      DATA ZILCH/0.0/, IONE/1/
C-----------------------------------------------------------------------
      INDX2(I,J,N)=I+(J-1)*N
C-----------------------------------------------------------------------
C     Set IFLAGS45 to specify type of guess.
C
C     IFLAG45 = 0 --- Read MOs from JOBARC
C     IFLAG45 = 1 --- Do a core guess
C     IFLAG45 = 2 --- NDDO guess
C     IFLAG45 = 3 --- Reserved (nothing implemented)
C     IFLAG45 = 4 --- Read MOs in SO basis from OLDMOS
C     IFLAG45 = 5 --- Read MOs in AO basis from OLDAOMOS
C     IFLAG45 = 6 --- Projected minimum basis set as guess
C     IFLAG45 = 7 --- Extended Huckel guess
C     IFLAG45 = 8 --- To override any defaults
C     IFLAG45 = 9 --- Atomic density superposition guess
C-----------------------------------------------------------------------
C
      ROHFMO     = .FALSE.
      G_FLE_READ = .FALSE.
      Geom_Opt   = (iFlags(h_IFLAGS_opt_method).GT.0)
      Findf_Vib  = (iFlags(h_IFLAGS_vib).GT.1) .or. 
     $     (iFlags2(h_IFLAGS2_prep_vibron).EQ.1) .or.
     $     (iFlags2(h_IFLAGS2_scan_normal).EQ.1)
      Brueckner  = (iFlags(h_IFLAGS_brueckner).GT.0)
      Rohf       = (iFlags(h_IFLAGS_reference).EQ.2)
      hfstab     = (iFlags(h_IFLAGS_hfstability).EQ.2)
      Occ_In_Zmt = (iFlags(h_IFLAGS_occupation).NE.0)
      IFLAG45    = iFLAGS(h_IFLAGS_guess)
      QRHF       = (iFlags(h_IFLAGS_qrhf_gen).GT.0)
      Length     = 0
      IKICK      = 0
      LUGSS      = 71
C
      IF      (IFLAG45 .EQ. 1) THEN
         GUESS_TYPE = "CORE"
      ELSE IF (IFLAG45 .EQ. 2) THEN
         GUESS_TYPE = "NDDO"
      ELSE IF (IFLAG45 .EQ. 3) THEN
         GUESS_TYPE = "EMPTY"
      ELSE IF (IFLAG45 .EQ. 4) THEN
         GUESS_TYPE = "EXPORT_MOS"
      ELSE IF (IFLAG45 .EQ. 5) THEN
         GUESS_TYPE = "EXPORT_AOS"
      ELSE IF (IFLAG45 .EQ. 6) THEN
         GUESS_TYPE = "PROJ_STO-3G"
      ELSE IF (IFLAG45 .EQ. 7) THEN
         GUESS_TYPE = "EHT" 
      ELSE IF (IFLAG45 .EQ. 8) THEN
         GUESS_TYPE = "OVERRIDE"
      ELSE IF (IFLAG45 .EQ. 9) THEN
         GUESS_TYPE = "ATOMIC"
      ELSE IF (IFLAG45 .EQ. 10) THEN
         GUESS_TYPE = "EXPORT_DENS"
      ELSE 
         GUESS_TYPE = "NONE"
      END IF

C Now, change the GUESS_TYPE SCF for runs that are not the very
C first one. Those include, Brueckner, hfsatbility, finite
C difference vibrational frequncy and geometry optimizations.
C
C For geometry optimizations we can read the previous step's
C vectors if they exsist (no change in symmetry from one step
C to the other).
C
      IF (Geom_Opt) THEN
         CALL GETREC(0,"JOBARC","SCFEVCA0",Length,iScr)
C
C If someone wants to always read from the guess file allow it.
C
         IF (GUESS_TYPE.EQ."EXPORT_MOS") THEN
            CALL READGS(NOCC,NIRREP,IUHF)
            IF (GSSALW) THEN
               GUESS_TYPE = "EXPORT_MOS"
C
C Allow read in vectros with out a GUESS file for geometry optimizations
C The occupations are from core Hamiltonian unless specified in ZMAT.

               IF (GSSOPT) G_FLE_READ = .TRUE.
            ELSE
               IF (Length.GT.0) GUESS_TYPE = "VECFRM_JARC"
            END IF
CSSS            G_FLE_READ = .TRUE.
         ELSE
            IF (Length.GT.0) GUESS_TYPE = "VECFRM_JARC"
         END IF

      END IF
C
C For finite difference calculations, lets allow the vectors 
C generated in the previous step to be used as the current steps 
C starting vectors via AOBASMOS file. The starting occupation 
C and the symmetry are all handled correctly. Another option is 
C to first genrate the starting vectors (OLDAOMOS) from the desired calculation
C and then use those as the starting vectors in subsequent finite
C difference calculation. The default is to use starting vectors
C from previous step. Ajith Perera, 07/07 

      AOBASMOS_FILE = .FALSE.
      AODENS_FILE   = .FALSE.
      INQUIRE(FILE='AOBASMOS',EXIST=AOBASMOS_FILE)
      INQUIRE(FILE='AODENS',EXIST=AODENS_FILE)
C
CSSS      AUTOMATE = GUESS_TYPE .NE. "OVERRIDE" .OR. 
CSSS     &           GUESS_TYPE .EQ. "NONE" 
C This is a good thing to do but symmetry handling is still somewhat
C at very preliminary stages and lots of jobs are failing. So, let's
C not use this for the moment. Ajith Perera, 05/08.
C
cmn
c This can be set on input now using the auto_occ keyword.
c
      if (IFLAGS2(h_IFLAGS2_auto_occ) .eq. 1) then
         AUTOMATE = .TRUE. 
      elseif (iFLAGS2(h_IFLAGS2_auto_occ) .eq. 0) then
         AUTOMATE = .FALSE.
      else
         write(6,*) ' AUTO_OCC flag in joda not processed '
         write(6,*) ' @initges: AUTOMATE not defined '
         call errex
      endif
#ifdef _DEBUG_LVLM1
      write(6,*) 'Initges: logicals '
      write(6,*) 'Findf_Vib : ', findf_vib
      write(6,*) ' aodens_file ', aodens_file
      write(6,*) ' aobasmos_file ', aobasmos_file
      write(6,*) ' Automate ', automate
#endif
c
      IF (Findf_Vib .AND. AODENS_FILE .AND. AUTOMATE) THEN
          IF (GUESS_TYPE .NE. "EXPORT_AOS") GUESS_TYPE = 
     &        "EXPORT_DENS"
      ELSE
          IF (GUESS_TYPE .EQ. "OVERRIDE") GUESS_TYPE = "CORE"
      ENDIF
C
C Brueckner case; in the near future the restriction on findif vibrational
C frq. is relaxed to incorporate the points for which the symmetry
C remain unchanged.
C
      BRUKMOS_EXIST = .FALSE.
      
Cmn      IF (Brueckner .and. .not.Findf_Vib) Then
      IF (Brueckner) Then
         inquire(file='BRUKMOS',exist=BRUKMOS_EXIST)
         iunit = 79
         CALL GETREC(0, "JOBARC", "SCFEVCA0", Length, Iscr)
         IF (Length.GT.0 .OR. BRUKMOS_EXIST) THEN
            GUESS_TYPE = "VECFRM_JARC"
         ELSE
            IF (GUESS_TYPE.EQ."NONE") GUESS_TYPE = "CORE"
         END IF
      END IF
C
C HFSTABILITY=FOLLOW case. I hope no body try to do Frequncies while
C following instabilities.
C
      IF (hfstab) Then

C Check to see whether this is a HFSTBILITY=FOLLOW calculation
C and actually following allowed instability (the second scf step).
C If that is the case, we should always read eigenvectos from
C the JOBARC file regardless of whether the other options are
C requested. Ajith Perera, 01/2001
C
         CALL GETREC(-1, "JOBARC", "SCFKICK ", IONE, IKICK)
         IF (IKICK.EQ.1) GUESS_TYPE = "VECFRM_JARC"
      END IF
C
C In principal the following call to READGS should be done by the
C the block that handle the GUESS and OLDMOS files. However, the
C GUESS file has a piece of info that is used to control whether
C the same starting OLDMOS file should be read during a subsequent
C SCF. Otherwise, the function of READGS is to read the GUESS file
C and initialize the common blocks /GSCOMA/ and /GSCOMB/.
C
C GSSOPT - .TRUE. if the GUESS file exists. If it is .TRUE. then the
C           initial guess options are read from the GUESS file.
C GSSALW - If it is .TRUE. then initial guess parameters are always
C          read from the GUESS file. If it is .FALSE. and GSSOPT is
C          .TRUE. then initial guess parameters are read from GUESS
C          for the first SCF calculation only.
C
      IF (.NOT.G_FLE_READ) THEN
          CALL READGS(NOCC,NIRREP,IUHF)
          IF (GSSALW) GUESS_TYPE = "EXPORT_MOS"
          IF (GSSOPT) G_FLE_READ = .TRUE.
C
      END IF
C
C If we are doing an ROHF calculation and there are no eigen-vectors
C in the JOBARC file, use the core Hamiltonian guess.
C
       IF (Rohf) THEN
          CALL ZERO(SCRN,NBAS*NBAS)
          CALL GETREC(-1,'JOBARC','SCFEVCAS',NBAS*NBAS*IINTFP,SCRN)
          ZTEST=DNRM2(NBAS*NBAS,SCRN,1)
C
          IF (ZTEST.GT.ZILCH) THEN
             IF (GUESS_TYPE.EQ."NONE") GUESS_TYPE = "VECFRM_JARC"
          ELSE
             IF (GUESS_TYPE.EQ."NONE") GUESS_TYPE = "CORE"
          END IF
       END IF

C
C If we have not choosen a procedure for initial until now, let's
C do it now. At the moment, we do the core guess. We can't do
C NDDO since NDDO execuation is controled by ME aces2 and we
C do not want to issue shell commands from all the other MEs. When
C the EH guess become reliable, we should change the core guess
C to EH. (*** old initges had checks for QRHF and vib***).
C
       IF (GUESS_TYPE.EQ."NONE") THEN

c       o default to core Hamiltonian
          GUESS_TYPE = "CORE"

c       o elevate to JOBARC records if possible
c          print *, '@INITGES: Initial occupation is'
c          print *, (NOCC(I),I=1,NIRREP),' /',(NOCC(J),J=9,8+NIRREP)
          J = 0
          DO I = 1, 16
             J = J + NOCC(I)
          END DO
c          CALL GETREC(0,"JOBARC","SCFEVCA0",Length,iScr)
c          IF (LENGTH.NE.0.AND.J.NE.0) GUESS_TYPE = "VECFRM_JARC"
          CALL GETREC(-1,"JOBARC","SCFEVCA0",nbas*nbas*iintfp,scr1)
          dtmp = ddot(nbas*nbas,scr1,1,scr1,1)
          IF (DTMP.GT.1.D-4.AND.J.NE.0) GUESS_TYPE = "VECFRM_JARC"
          
       END IF
C
C Get the starting occupation from the core Hamiltonian unless
C
C      1. The user has specified it in the ZMAT
C      2. The user has specified it in the GUESS file to go along
C         with OLDMOS file
C      3. Provided the AOBASMOS file. In which case the occupation
C         is in the AOBASMOS file.
C
      IF (.NOT. (Occ_In_Zmt
     &    .OR.  (GUESS_TYPE .EQ. "EXPORT_MOS ")
     &    .OR.  (GUESS_TYPE .EQ. "EXPORT_AOS ")
     &    .OR.  (GUESS_TYPE .EQ. "VECFRM_JARC" ))) THEN
C
         CALL GEN_CORE_HAM_OCC(SCR1, SCRN, ONEH, EVAL, SCRA,
     &                         SCR2, SCRB, XFORM, SCRTMP,
     &                         NOCC, LDIM1, LDIM2, IUHF, NBAS)
C
C  If the user DID specify a guess, then lock this occupation.
C
CSSSS      ELSE IF (GUESS_TYPE.EQ."VECFRM_JARC".AND..NOT.Occ_In_Zmt) THEN
      ELSE IF (GUESS_TYPE.EQ."VECFRM_JARC") THEN
C
C If we are here, we know that we have eigenvector JOBARC record. Then
C we must know something about the occupations as well. The occupation
C is read in to NOCC alredy by the main program. So, we do not
C need to do anything here.
C
C Ammendment: If QRHF is being used along with Breuckner orbital
C option set, what is read in the main program correspond to the
C post QRHF occupation. We need pre-QRHF occupation at this point.
C These occupation numbers are saved in QRHFIT. A. Perera, 10/2013.
C At this point, we have a occupation (either user specified it
C in the zmat or generated in the previous cycle, so we need to
C lock it. Otherwise tests like zmat.013b will fail. A. Perera, 12/2013. 
 
           IF (Brueckner .AND. QRHF) THEN
              CALL GETREC(20, "JOBARC", "ORGOCCA", 8, NOCC(1))
              CALL GETREC(20, "JOBARC", "ORGOCCB", 8, NOCC(9))
           ENDIF
                              
           IFLAGS(h_IFLAGS_lock_orbocc) = 1

CSSS         CALL GETREC(20,'JOBARC','OCCUPYA0',NIRREP,NOCC(1))
CSSS         IF (IUHF.EQ.1) THEN
CSSS            CALL GETREC(20,'JOBARC','OCCUPYB0',NIRREP,NOCC(9))
CSSS         ELSE
CSSS            CALL ICOPY(8,NOCC(1),1,NOCC(9),1)
CSSS         END IF
CSSS         WRITE(6,5001)(NOCC(I),I=1,NIRREP)
CSSS         WRITE(6,5002)(NOCC(8+I),I=1,NIRREP)
CSSS         WRITE(6,5010)
C
      ELSE
C
C If we are here we need to check whether we have a occupation.
C If the occupation is read from AOBASMOS or ZMAT, we need to
C lock it. If the OLDMOS file read along with the GUESS file
C then we also have a occupation. But, if the user choose GUESS
C READ_SO_MOS but did not file a GUESS file, we should be
C able to continue, but needs an occupation. Let's generate it by
C using the core Hamiltonian.
C
C The final occupation numbers were added to the NEWMOS/OLDMOS file.
C So, we no longet needed to generate core Hamiltonian occupation 
C for EXPORT_MOS option. Ajith Perera 10/2015.
C
CSSS         IF (GUESS_TYPE.EQ."EXPORT_MOS".AND..NOT.Occ_In_Zmt.AND.
CSSS     &       .NOT.G_FLE_READ) THEN
C
         IF (.NOT.Occ_In_Zmt.AND. .NOT.G_FLE_READ .AND. 
     &       .NOT. (GUESS_TYPE.EQ."EXPORT_MOS")) THEN
            CALL GEN_CORE_HAM_OCC(SCR1,SCRN,ONEH,EVAL,SCRA,SCR2,SCRB,
     &                            XFORM,SCRTMP,NOCC,LDIM1,LDIM2,IUHF,
     &                            NBAS)
         END IF
C
C We have not get any info from the GUESS file (it is not there).
C So we must set the GSCOMA common block variables. However, note that
C only thing we can alow is just to read the orbitals; no locking, swaping
C rhf guess from uhf, etc. If you want them, you have to put a
C GUESS file.
C
         GSSRED = .TRUE.
         IF (Occ_In_Zmt) IFLAGS(h_IFLAGS_lock_orbocc) = 1
C
      END IF
C
      CALL ICOPY(16,NOCC,1,IOLDOCC,1)
C
C The core Hamiltionian guess. During a geometry optimization
C use the last step's vectors from the "JOBARC" file (ie. core guess
C is done only in the  zero geo, opt. cycle. The same can
C be done for findif vibrational frequncy calculations as long
C as the symmetry is unchanged between the current and the previous point.
C At the moment this is not done since the JOBARC truncation
C wipes out the previous vectors. As soon as we are convinced that
C the code functions properly when trunction is done only when the
C point group is changed, we are able to use the eigenvectros
C of the previous step stored in the JOBARC file. In general
C it is possible to use vectros from the previous step regardless
C of the symmetry as long as the appropriate transformations are
C done (a crude attempt at this is used in the option
c GUESS_TYPE=EXPORT_AOS, see below). The HFSTABILITY=FOLLOW
C must also use the exsisting eigenvectors from the JOBARC.
C
      IF (GUESS_TYPE.EQ."CORE") THEN
C
         CALL CORE_HAMILTON_GUESS(SCR1,SCRN,ONEH,SCRA,SCR2,SCRB,XFORM,
     &                            SCRTMP,LDIM1,LDIM2,IUHF,NBAS)
C
C The NDDO guess. All the preceeding comments for core guess also
C applies to the NDDO guess. One important difference though is that the
C NDDO guess is generated by a seperate ME that run prior to
C the vscf ME. It is controled by the driver xaces2 (aces2).
C For geo. opts. the nddo runs only in the zero cycle. For the finite
C difference Freq. at the moment, every finite difference point
C needs to run the NDDO. It will soon be changed to only when the
C comp. point group change. The HFSTABILITY=FOLLOW must use the
C exsisting vectors.
C
      ELSE IF (GUESS_TYPE.EQ."NDDO") THEN
C
         CALL NDDO_GUESS(DENS, LDIM1, IUHF)
C 
C The neutral atomic density superposition guess for the density
C uses the previously generated spin and spatially averaged atomic densities
C for the ground state multiplicity of each atom in the molecule and
C generates a block diagonal density matrix.  This matrix is not idempotent.
C         
      ELSE IF (GUESS_TYPE.EQ."ATOMIC") THEN
C       
         CALL ATOMIC_GUESS(DENS, LDIM1, IUHF)
C          
C The extended Huckel or projected minimal basis guess. All the
C preceeding comments for core guess applies. The EHT and
C projected STO-3G basis guess is generated internally.
C For geo. opts. these guesses generated only in the zero cycle.
C For the finite difference Freq. at the moment, every finite difference
C point needs to generate EHF/P-STO-EG guess. It will soon be changed
C to only when the comp. point group change. The HFSTABILITY=FOLLOW
C must use the exsisting vectors.

      ELSE IF (GUESS_TYPE.EQ."EHT".OR.GUESS_TYPE.EQ."PROJ_STO-3G") THEN
C
         CALL EHT_PRJSTO3G_GUESS(DENS,CZAO,CSO,ZAOSO,SCR1,SCRN,DCORE,
     &                           IUHF,NBAS,NBASX,LDIM1)
C
C The GUESS_TYPE ="EXPORT_AOS  " applies only for finite difference
C vibrational frequency calculations. It is simply an extension of
C the EXPORT_MOS option. During finite difference frequency calcuations
C the computational point group changes (lowered) and the straight forward
C use of SO basis eigenvectors in the JOBARC or OLDMOS from previous
C point is not possible unless the symmetry remain unchanged. The
C EXPORT_AOS use the MO (converted to AO basis) file generated from a
C separte SCF calculation and it is being read for every finite difference
C step (the same vectors)(NO NDDO, EHT or any other guess when this option
C is choosen).
C
      ELSE IF (GUESS_TYPE.EQ."EXPORT_AOS") THEN
C
C These tests should be performed in JODA. Pepole should not come
C this far to figure out that they can not use this option with
C out GUESS file and OLDMOS file. I leave them here only as a duplication.
C *****Anthony please trap in joda.
C
         OLDAOMOS_FILE = .FALSE.
         INQUIRE(FILE='OLDAOMOS',EXIST=OLDAOMOS_FILE)
         IF (AOBASMOS_FILE .OR. OLDAOMOS_FILE) THEN
C
C Check and see whether we are following a instability. These are
C part of the old versions's checks, and in principle no longer
C needed. I will leave them just in case some one find them useful.

            IF (IKICK .EQ. 1) THEN
               WRITE(6,*) "@INITGES: GUESS=READ_AO_MOS and ",
     &                    "HFSTABILITY=FOLLOW are incompatible."
               CALL ERREX
            END IF
         ELSE
            WRITE(6,*) '@INITGES: File OLDAOMOS/AOBASMOS needed for',
     &                 ' GUESS=READ_AO_MOS'
            CALL ERREX
         END IF
C
         CALL EXPORT_GUESS_AO_ORBS(CZAO, EVAL, OCCNUM, REFCOORD,
     &                             COORD, CSO, ANGTYP, ATMASS,
     &                             CNEW, DCORE, ZAOSO, SCR1, SCRN,
     &                             SCRA, SCRB, SCRTMP, IUHF, NBAS,
     &                             NBASX, MAXDCOR, LDIM2, NATOMS,
     &                             IOLDOCC)
C
C The GUESS_TYPE = "EXPORT_MOS" import symmetry adapated eiegnvectors
C from a previous calculation. The file "GUESS" is needed to manipulate
C how the guess vectors are being used. During a finite difference
C frequency calculation or geometry optimization (analytical or numerical)
C the "GUESS" file and the eigenvectors are read only once in the
C first cycle. In the subsequent cycles we can use the vectors from
C the previous cylces only for geometry optimizations (NO NDDO or EHT
C or any other guess when this option is choosen).
C
      ELSE IF (GUESS_TYPE .EQ. "EXPORT_MOS ") THEN
C
C These tests should be performed in JODA. Pepole should not come
C this far to figure out that they can not use this option with
C out GUESS file and OLDMOS file. I leave them here only as a duplication.
C ******Anthony please trap in joda.
C The GSSOPT get set in READGS and set to .TRUE. when the GUESS file
C is present.
C
         IF (.NOT. GSSOPT .AND. GSSALW)THEN
            WRITE(6,*) '@INITGES: File GUESS needed for ',
     &                 'GUESS=READ_SO_MOS'
            CALL ERREX
         END IF
C
         OLDMOS_FILE = .FALSE.
         INQUIRE(FILE='OLDMOS',EXIST=OLDMOS_FILE)
C
         IF (.NOT. OLDMOS_FILE)THEN
            WRITE(6,*) '@INITGES: File OLDMOS needed for ',
     &                 'GUESS=READ_SO_MOS'
            CALL ERREX
         END IF

         CALL EXPORT_GUESS_MO_ORBS(SCR1,SCRN,SCRA,SCRB,SCR2,LDIM2,NBAS,
     &                             IUHF,ROHF)

      ELSE IF (GUESS_TYPE .EQ. "VECFRM_JARC") THEN
C
C Read from JOBARC when we are in the second macro-iteration of a
C Brueckner iteration. This is signified by IFLAGS(16)=0 (i.e. 0
C SCF cycles). This SHOULD work in finite difference jobs too.
C In the first macro-iteration of a finite difference Brueckner
C job, a standard guess is mades (core, nddo, ....).
C
         IF (BRUKMOS_EXIST) THEN
C
C The BRUKMOS file is created in mrcc. This is one of Marcel
C things and I fail to see what he trying to acomplish. If all
C he wants to read from the previous vectors in JOBARC it is
C done already witout BRUKMOS.
C
            perline = 4
            call get_matrix(scr1, nbas, nbas, perline, iunit)
            call putrec(20,'JOBARC','SCFEVCA0',nbas*nbas*iintfp, scr1)
            call putrec(20,'JOBARC','SCFEVCB0',nbas*nbas*iintfp, scr1)
C
C IFlags(16): SCF_MAXCYC  = 0
C IFlags(38): NON_HF      = 1
C IFlags(34): QRHF_ORBITAL= 0
C IFlags(77): QRHF_GENERAL= 0
C
            iflags(16)=0
            iflags(38)=1
            iflags(34)=0
            iflags(77)=0
            call putrec(20, 'JOBARC', 'IFLAGS  ',100,iflags)
c
c also construct reordering vectors as in vtran
c
            CALL GETREC(20,'JOBARC','OCCUPYA ',NIRREP, NOCC_IRREP)
C
            ICNT=0
            CALL IZERO(IMAGE,NBAS)
            DO I=1,NIRREP
               DO J=1,NBFIRR(I)
                  ICNT=ICNT+1
                  IMAGE(ICNT)=I
               END DO
            END DO
C
            CALL FMREVC(scrn,NOCC_IRREP,IMAGE,nbas,NIRREP,1,IUHF,"VSCF")
C
          END IF
C
          CALL EIGVECTORS_FRM_JOBARC(SCR1,SCRN,SCRA,SCRB,SCRTMP,LDIM2,
     &                               NBAS,IUHF)
C
C QRHF or TD-CCSD Brueckner requires undoing the swapping of orbital
C that happens after the scf (during creating new reference using QRHF
C QRHF_GEN,QRHF_ORB,QRHF_SPIN keywords), 11/2105, Ajith Perera.
C
          If (iFlags(h_IFLAGS_qrhf_gen).NE.0)
     &       CALL FIX_QRHF_EVECS(SCR1,SCRN,DCORE,MAXDCOR,NBAS,IUHF)

      ELSE IF (GUESS_TYPE .EQ. "EXPORT_DENS") THEN
          
          CALL EXPORT_AODENS(DCORE, DENS, SCR2, LDIM1, LDIM2, MAXDCOR, 
     &                       NBASX, NBAS, IUHF)
C
      ELSE  

         IF (IFLAG45.EQ.3)THEN
            WRITE(6,*) '@INITGES: GUESS option is not available.'
            CALL ERREX
         END IF
C
      END IF
C
C-----------------------------------------------------------------------
C
C Make the density matrix for all the methods except for NDDO,
C EHT and projected STO-3G guess. For those cases the density
C has already been made. Also, for ROHF cases also makes the
C contribution to the full guess eigenvectors from this sym. blcok.

      CALL MAKE_GUESS_DEN(SCR1, SCR2, SCRA, SCRN, DENS, NBAS, LDIM1,
     &                    LDIM2, IUHF, ROHF, ROHFMO, GUESS_TYPE)

      CALL GETREC(-1,'JOBARC','OSCALC  ',1,IOS)
      If (IOS .gt. 0 .AND. GUESS_TYPE .EQ. "VECFRM_JARC") 
     &    CALL UNO_4TDCC(DENS,XFORM,EVAL,CSO,DCORE,MAXDCOR,LDIM1,
     &                    ISQRLN(NIRREP+1),IUHF,NBAS)

      RETURN
      END

