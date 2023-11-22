
      PROGRAM VSCF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER(LUINT=10)
      PARAMETER(LUNITSE=25)
      PARAMETER(LUHF2=15)
#include "aces.h"
#include "flags.h"
#include "kscore.com"
#include "molcas.com"
      INTEGER MatDim(8,2)

       logical scfksiter,kshf,scfks,scfksexact,scfkslastiter,agrad
       logical hfdftgrad, hfdft,oepKS, COT, DENS, ORBS, props
       integer setptr
c
#include "iuhf.com"
#include "intgrtflags.com"
#include "oep.com"

      LOGICAL transfer
      logical slater,loslater
      CHARACTER*3 LRPP(3)
      CHARACTER*4 CREF(3),CREFJ,FULGRP,CMPGRP
      CHARACTER*8 TITLE(24), LABEL
      CHARACTER*23 CGUESS(8)
      CHARACTER*5 PROGRAM 
      CHARACTER*80 FNAME
      INTEGER SCFPRNT
      INTEGER SWAP,READMO,WRITMO
C
C For HF-DFT Gradients, P. Verma, A. Perera 08/2008 
C
      integer ioccp(2),ivrtp(2),ioccvrt(2)
C
c...watson...
c   BEGIN: Putting in the necessary variables for the DKH routine.
      logical dkhess
      integer dkhorder,DKHZMAX
c...watson...

      LOGICAL AOFIL, POSTSCF
      LOGICAL GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      LOGICAL ROHFMO
      LOGICAL NOCONV,YESNO,YESNO2,CHGOCC,ISTHERA,ISTHERB,TERMINATE
      logical dmpflg
      LOGICAL QC_PHASE,BRUECKNER
      LOGICAL CHKPHASE
      LOGICAL FIRSTSCF,YESNO3, bTmp, IIII_EXIT
      DOUBLE PRECISION MAX_ERROR
      DIMENSION IOLDOCC(16)
      DIMENSION IRPSZA(37),IRPSZB(30)
      DIMENSION NOCC(16)
      DIMENSION NDET(16)
      DIMENSION NPRVOCC(16)
      DIMENSION IADERR(200),IADFCK(200,2)
C
      COMMON // ICORE(1)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/  LUOUT,MOINTS
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /POPUL/  NOCC
      COMMON /FOCK/   AOFIL
#include "symm2.com"
      COMMON /RPPINF/ ISTART,IORDER,IADERR,IADFCK
      COMMON /KILLS/  KILLFLG(3)
      COMMON /SCFPRINT/ SCFPRNT
      COMMON /GSCOMA/ GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      COMMON /GSCOMB/ SWAP(4,8,2),LOCK(8,2),IPRTGS(8,2),ISTOP1,ISTOP2,
     1                READMO,WRITMO,IUHFRHF,LUGSS
      COMMON /GSCOMC/ ROHFMO
      COMMON /DMPCOM/ ITDLST
      COMMON /LSHCOM/ ALPHA1,BETA1
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW

      COMMON /LOGBLOCK/ ITER

#include "molcas.dat"
cwc1
#include "mpif.h"
#include "parallel_aces.com"
      logical bExist,dirscf,fdiff,schwrz,pople
      common /GMSdirect/ dirscf,fdiff,schwrz,pople,i011,i012,i013
      character*80 str
      data one,zer /1.d0,0.d0/
cwc0      

      DATA CREF /' RHF',' UHF','ROHF'/
      DATA CGUESS /'MOREAD                 ',
     &             'CORE                   ',
     &             'NDDO                   ',
     &             'WALT_PRJDEN (NCI)      ',
     &             'READ_SO_MOS            ',
     &             'READ_AO_MOS            ',
     &             'MIN_BASIS              ',
     &             'HUCKEL                 '/
      DATA LRPP /'OFF','RPP',' QC'/
      DATA IONE /1/
      DATA I16 /16/

c   o variables for QC-SCF
      DATA NDET /16*0/
      DATA NPRVOCC,IPREVOCC /16*0,0/
      DATA QC_PHASE /.FALSE./
      DATA QC_THRESH /1.d-2/

         call callstack_init('VSCF')
          
c ----------------------------------------------------------------------

c   o initialize the MPI environment
      call MPW_Init(iErr)
      if (iErr.ne.MPI_SUCCESS) then
         print *, '@VSCF: unable to initialize the MPI environment'
         call aces_exit(1)
      end if

c   o initialize the parallel_aces common block
      call aces_com_parallel_aces

#ifdef _HAVE_MPI
c   o change to the temp directory
      call cd_noderank(irank,iErr)
      if (iErr.ne.0) then
         print *, '@VSCF: unable to cd to the temp directory'
         print *, '       attempting to run in current dir'
      end if
#endif

c ----------------------------------------------------------------------

      LUOUT=6
      MOINTS=50
      YESNO=.FALSE.
      YESNO2=.FALSE.
      TERMINATE=.FALSE.

      CALL IZERO(KILLFLG,3)
      CALL TIMER(1)
      CALL CRAPSI(ICORE(1),IUHF,-1)
cwc1
      dirscf=.false.
      fdiff=.false.
      schwrz=.false.
      pople=.false.
      scfks=.false.
      kshf =.false.
      scfksexact=.false.
      hfdftgrad=.false.
      agrad=.false.
      hfdft=.false.
      loslater=.false.

      IF (iflags(h_IFLAGS_integrals) .EQ. 5) THEN
      IF (iflags2(h_IFLAGS2_direct) .EQ. 1) THEN
         dirscf=.true.
      ENDIF
      ENDIF
      IF (iflags2(h_IFLAGS2_program) .EQ. 3) THEN
         PROGRAM = "MN_A3"
      ENDIF 
c
c-----initilize the oep timers------------
      if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4) then
         call oeptimer ('', 'INI')
         call oeptimer ('VSCF', 'ON')
         call oeptimer ('VSCF: PRE-SCF', 'ON')
      endif
c 
      if (dirscf) then
         open (98,file='ZMAT',status='OLD')
         str=' '
         i=1
         do while (str.ne.'*GAMESS'.and.i.ne.0)
            read (98,'(a)',end=8) str
            goto 9
 8          i=0
 9          continue
         end do
         if (i.ne.0) read (98,'(a)') str
         do while (str.ne.' '.and.str(1:1).ne.'*'.and.i.ne.0)
            i=1
            do while (str(i:i).eq.' '.and.i.lt.(len(str)-8))
               i=i+1
            end do
            if (str(i:i+6).eq.'SCHWRZ=') then
               if (str(i+7:i+7).eq.'T') schwrz=.true.
            else if (str(i:i+5).eq.'FDIFF=') then
               if (str(i+6:i+6).eq.'T') fdiff=.true.
            else if (str(i:i+5).eq.'POPLE=') then
               if (str(i+6:i+6).eq.'T') pople=.true.
            else if (str(i:i+6).eq.'MEMORY=') then
c              ! compatibility w/ integral derivatives
            else if (str.ne.' '.and.str(1:1).ne.'*') then
               if (irank.eq.0) then
                  write (*,*) 'error in *GAMESS namelist'
                  write (*,*) 'Each line should be of the form:'
                  write (*,*) '(SCHWRZ|FDIFF|POPLE)=(T|F)'
               end if
               call aces_exit(1)
            endif
            read (98,'(a)',end=10) str
            goto 11
 10         i=0
 11         continue
         end do
         close (98)
      end if

      if (.not.dirscf) schwrz=.false.
      if (.not.schwrz) fdiff=.false.
c temporary workaround
      pople=.false.

      IF (dirscf) then
         open (98,file='MOL')
         if (irank.eq.0) open (99,file='GAMESS.LOG')
         call GETINF (98,.true.,99)
         close (98)
         if (pople) then
            intg76=1
         else
            intg76=0
         end if
         if (irank.eq.0) then
            write (*,*) 'GAMESS initialization complete'
            write (99,'(a)') 
     >     'direct SCF mode, two-electron GAMESS integrals will be used'
            if (schwrz) then
               write (99,'(a)') 'Schwarz estimator is ON'
            else
               write (99,'(a)') 'Schwarz estimator is OFF'
            end if   
            if (fdiff) then
               write (99,'(a)') 'incremental Fock matrix will be built'
            else
               write (99,'(a)') 'full Fock matrix will be built'
            end if   
         end if   
      END IF
cwc0
      CALL IZERO(NOCC,16)
      SCFPRNT=IFLAGS(1)
      MAXIT=IFLAGS(16)
      IEXP=-IFLAGS(5)
      TOL=10.**IEXP
      MAXCOR=IFLAGS(36)
      IBUFLN=600
C
CJDW 3/28/97. Check whether this is the first entry to xvscf.
C             FIRSTSCF is used to control whether we print certain
C             things or not.
C
      CALL GETREC(1,'JOBARC','FIRSTRUN',1,ITMP)
      FIRSTSCF = (itmp.ne.0).or.(iFlags(1).ge.10)

cmn put flag on JOBARC for later use
      CALL PUTREC(20,'JOBARC', 'STEOM_SZ', IONE, 0)
C
C AOFIL IS THE LOGICAL VARIABLE WHICH DECIDES WHETHER FOCK MATRICES
C ARE CONSTRUCTED FROM PK LIST OR AO INTEGRALS
      AOFIL=(IFLAGS(96).NE.0)

C     ----- Set Damping parameters -----
      IF (IFLAGS2(109).EQ.1) THEN
         DMPFLG = .TRUE.
         DMPTOL = 0.01D+00 * DFLOAT(IFLAGS2(110))
         DAMP   = DMPTOL + 0.01D+00
         DEAVG  = 0.0D+00
      ELSE IF (IFLAGS2(109).EQ.2) THEN
         DMPFLG = .TRUE.
         DMPTOL = 0.01D+00 * DFLOAT(IFLAGS2(110))
         DAMP   = 0.01D+00 * DFLOAT(IFLAGS(15))
         IF (DAMP.LT.DMPTOL) THEN
            write(6,*) ' @VSCF-F, Static damping parameter is < DMPTOL '
            CALL ERREX
         END IF
      ELSE
         DMPFLG = .FALSE.
         ITDLST = 0
      END IF

C     ----- Set level-shifting parameters -----
      IF (IFLAGS2(111).EQ.0) THEN
         ALPHA1 = 0.0D+00
      ELSE
         ALPHA1 = 0.01D+00 * DFLOAT(IFLAGS2(111))
      END IF
      IF (IFLAGS2(111).EQ.0) THEN
         BETA1 = 0.0D+00
      ELSE
         BETA1 = 0.01D+00 * DFLOAT(IFLAGS2(112))
      END IF

C  LET'S GET SOME SYMMETRY INFO FROM THE INTEGRAL FILE SO WE CAN
C  CALCULATE SIZES OF ARRAYS AND WHAT NOT.

      SEWARD=(IFLAGS(56).EQ.4)

      if (seward) then

C Get the information about the format of the two-electron integrals file. This 
C only need when we use seward to generate integrals. SEWARD can store 
C integrals in old VMOL format (IIII, IJIJ ..) or it can store integrals
C in petite list format. We can use both storage schemes. The new routines
C were added to handle the petite list format (Anthony Yau Aug. 2000).
         INQUIRE(FILE='ORDINT',EXIST=PETITE_LIST)

         call gfname('ONE_INT ',fname,ilength)
         inquire(file=fname(1:ilength),exist=yesno2)
         if (yesno2) then
            open(lunitse,file=fname(1:ilength),
     &           form='UNFORMATTED',access='SEQUENTIAL')
            read(lunitse) title, nirrep,
     &                    (nbfirr(i),i=1,nirrep),
     &                    repuls, iijunk
            close(lunitse,status='KEEP')

            if (petite_list) then
c            o start the MOLCAS backend
cSEWARD               call start("VSCF")
cSEWTRACE               call qenter("VSCF")
#ifdef _TRAPFPE
               call ignorefpe
cSEWARD               call sysinf
               call trapfpe
#else
cSEWARD               call sysinf
#endif /* _TRAPFPE */

c            o open the 2-electron integral file
               itmp=-1
               iopt=0
               call opnord(itmp,iopt,fnord,luord)
               if (itmp.ne.0) then
                  write(*,*) "@VSCF: ERROR - The MOLCAS 2-electron ",
     &                       "integral file could not be opened."
                  call errex
               end if

c            o initialize array pointers for Seward integral subroutines
c              (cf. $MOLCAS/src/scf/setup.f)...
               ipmat(1,lbbt) = 1
               do itmp = 1,(nirrep-1)
                  ipmat(itmp+1,lbbt)=ipmat(itmp,lbbt)
     &                               +nbfirr(itmp)*(nbfirr(itmp)+1)/2
               end do
               ipmat(1,lbbs) = 1
               do itmp = 1,(nirrep-1)
                  ipmat(itmp+1,lbbs)=ipmat(itmp,lbbs)
     &                               +nbfirr(itmp)*nbfirr(itmp)
               end do

c           end if (petite_list)
            end if

c        else if [ONE_INT does not exist]
         else
            write(*,*) "@VSCF: Unable to open an AO integral file!"
            call errex
         end if

c     else if (.not.seward) then
      else
         PETITE_LIST = .FALSE.
         CALL GFNAME('IIII    ',FNAME,ILENGTH)
         INQUIRE(FILE=FNAME(1:ILENGTH),EXIST=YESNO2)
         IF (YESNO2) THEN
            OPEN(LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &           ACCESS='SEQUENTIAL')
            READ(LUINT)TITLE,NIRREP,(NBFIRR(I),I=1,NIRREP),REPULS,
     &                 IIJUNK
            CLOSE(LUINT,STATUS='KEEP')
         ELSE
            WRITE(LUOUT,*)
     &         "@VSCF-F, Unable to open an AO integral file!"
            CALL ERREX
         END IF

c     end if (seward)
      end if

C Dump out the nuclear repulsion energy from the integral program. This
C is the correct energy to use in any subsequent calculations, since it
C is the one which goes with the integrals and the SCF energy.
      CALL PUTREC(20,'JOBARC','NUCREP  ',IINTFP,REPULS)
      NBAS=0
      DO I=1,NIRREP
         NBAS=NBAS+NBFIRR(I)
      END DO
C
      IF (FIRSTSCF .OR. IFLAGS(1).GE.1) THEN
         WRITE(LUOUT,'(a,i4,a,/)') "@VSCF: There are ",
     &      NBAS," functions in the AO basis."
         WRITE(LUOUT,'(a,i4,a,/)') "       There are ",
     &      NIRREP," irreducible representations."
         WRITE(LUOUT,'(/,T8,"Irrep",8X,"# of functions")')
         DO I=1,NIRREP
            WRITE(LUOUT,'(T10,I1,16X,I3)')I,NBFIRR(I)
         END DO
         WRITE(LUOUT,*)
      END IF
C
      IF (IFLAGS(77).NE.0) THEN
         CREFJ="QRHF"
      ELSE
         CREFJ=CREF(IFLAGS(11)+1)
      END IF
C
C  Handle all the logic for getting the correct initial guess,
C  whether it be reading the MO's and getting the occupation from
C  them, or using them as a starting point, etc...
C  If the ORGOCC (A or B) records exist, then the JOBARC correspond
C  to QRHF calcualtion. The occuptaion must read from those records.
C  09/2015, Ajith Perera.
C 
      CALL GETREC(0,'JOBARC','ORGOCCA',Length,NOCC(1))

      IF  (Length .GT. 0) THEN
          CALL GETREC(20,'JOBARC','ORGOCCA',NIRREP, NOCC(1))
      ELSE
          CALL GETREC(-1,'JOBARC','OCCUPYA0',NIRREP,NOCC(1))
      ENDIF 
      
      IF (IUHF.EQ.1) THEN
         CALL GETREC(0,'JOBARC','ORGOCCB',Length,NOCC(9))
         IF (Length .GT. 0) THEN 
             CALL GETREC(20,'JOBARC','ORGOCCB',NIRREP, NOCC(9))
         ELSE
             CALL GETREC(-1,'JOBARC','OCCUPYB0',NIRREP,NOCC(9))
         ENDIF 
      ELSE
         CALL ICOPY(8,NOCC(1),1,NOCC(9),1)
      END IF

      ISTHERA=.FALSE.
      ISTHERB=.FALSE.
      DO I=1,8
         IF (NOCC(I).NE.0) ISTHERA=.TRUE.
         IF (NOCC(8+I).NE.0) ISTHERB=.TRUE.
      END DO
      IF (.NOT.ISTHERA.OR.(IUHF.EQ.1.AND..NOT.ISTHERB).OR.
     &    IFLAGS(45).EQ.2) THEN
         CALL GETREC(20,'JOBARC','OCCUPYA ',NIRREP,NOCC(1))
         IF (IUHF.EQ.1) THEN
            CALL GETREC(20,'JOBARC','OCCUPYB ',NIRREP,NOCC(9))
         ELSE
            CALL ICOPY(8,NOCC(1),1,NOCC(9),1)
         END IF
         ISTHERA=.FALSE.
         ISTHERB=.FALSE.
         DO I=1,8
            IF (NOCC(I).NE.0) ISTHERA=.TRUE.
            IF (NOCC(8+I).NE.0) ISTHERB=.TRUE.
         END DO
      END IF
C
      CALL GETCREC(20,'JOBARC','FULLPTGP',4,FULGRP)
C
C  Check to see if this molecule is linear or not.
C
      IF (FULGRP(1:3).EQ.'D8h'.OR.FULGRP(1:3).EQ.'C8v') THEN
         CALL GETREC(20,'JOBARC','LINEAR  ',IONE,ILIN)
         IF (ILIN.EQ.1) FULGRP(2:2)='X'
      END IF
C
      CALL GETCREC(20,'JOBARC','COMPPTGP',4,CMPGRP)
      IF (FIRSTSCF .OR. IFLAGS(1).GE.1) THEN
         WRITE(LUOUT,5301) CREFJ,IFLAGS(16),FULGRP,CMPGRP
         IF (IFLAGS(45).EQ.0.AND..NOT.ISTHERA.AND.(IUHF.EQ.0.OR.
     &       .NOT.ISTHERB)) THEN
            WRITE(LUOUT,5307)CGUESS(2)
         ELSE
            WRITE(LUOUT,5307)CGUESS(IFLAGS(45)+1)
         END IF
         WRITE(LUOUT,5304) IFLAGS(5),LRPP(IFLAGS(10)+1),IFLAGS(40),
     &                     IFLAGS(20)
      END IF
C
      IF ((IFLAGS(17).NE.0.OR.(ISTHERA.AND.ISTHERB)).AND.
     &    IFLAGS(45).NE.1) THEN
         WRITE(LUOUT,5302)(NOCC(I),I=1,NIRREP)
         WRITE(LUOUT,5303)(NOCC(8+I),I=1,NIRREP)
         WRITE(LUOUT,5305)
      ELSE IF (IFLAGS(45).EQ.1) THEN
         WRITE(LUOUT,5306)
      END IF
C
      IF (FIRSTSCF .OR. IFLAGS(1).GE.1) THEN
         IF (IFLAGS2(109).GT.0) THEN
            IF (IFLAGS2(109).EQ.1) THEN
               WRITE(LUOUT,5320) DMPTOL
            ELSE
               WRITE(LUOUT,5330) DAMP
            END IF
         END IF
         IF (IFLAGS(11).EQ.2) THEN
            WRITE(LUOUT,5340) ALPHA1,BETA1
         END IF
      END IF

 5301 FORMAT(T3,'@VSCF-I, Parameters for SCF calculation: ',/,
     &       T8,'      SCF reference function: ',A4,/,
     &       T8,'Maximum number of iterations: ',I4,/,
     &       T8,'   Full symmetry point group: ',A4,/,
     &       T8,'   Computational point group: ',A4)
 5307 FORMAT(T8,'      Initial density matrix: ',A23)
 5302 FORMAT(T8,'   Alpha population by irrep: ',8(I3,2X))
 5303 FORMAT(T8,'    Beta population by irrep: ',8(I3,2X))
 5305 FORMAT(/)
 5306 FORMAT(T8,' *** Initial occupancies estimated from core ',
     &          'Hamiltonian. ***',/)
 5304 FORMAT(T8,'   SCF convergence tolerance: 10**(-',I2,')',/,
     &       T8,'    Convergence acceleration: ',A3,/,
     &       T8,'        Latest start for RPP: ',I3,/,
     &       T8,'                   RPP order: ',I3,/)
 5320 FORMAT(T8,'  Davidson dynamical damping: ON ',/,
     &       T8,'  Damping tolerance         ',F5.2)
 5330 FORMAT(T8,'              Static damping: ON ',/,
     &       T8,'              Damp factor   ',F5.2,/)
 5340 FORMAT(T8,' Doubles-singles  level-shifter : ',F6.2,/,
     &       T8,' Singles-virtuals level-shifter : ',F6.2,/)
C
      IF (IFLAGS(5).LT.7) THEN
         WRITE(LUOUT,*) "@VSCF-W, SCF convergence tolerance of ",
     &      "10**(-",IFLAGS(5),") is not too stringent."
         WRITE(LUOUT,
     &      '(T12,"*****   Use results with CAUTION!!!   *****")')
      END IF
C
C Compute the size of the single excitation matrix for QC-SCF. These
C sizes are exact only when there is a user-specified occupation. Since
C the occupation can change during the iterations and since it is safer
C to over-allocate than to under-allocate, we assume that the occupation
C per irrep can be changed by +/- 3. In the case of no occupation at
C all, let's allocate what is needed for the non-symmetry case using the
C total number of electrons and the basis functions.
C Ajith Perera and Carlos Taylor 04/06.
C
      NDIMRPA = 0
      IF (IFLAGS(h_IFLAGS_scf_extrap).EQ.2) THEN
         IF (ISTHERA.OR.ISTHERB) THEN
            DO IRREP = 1, NIRREP
               NOCCA = NOCC(IRREP)
               NOCCB = NOCC(8+IRREP)
               NVIRTA  = NBFIRR(IRREP)-NOCCA
               NVIRTB  = NBFIRR(IRREP)-NOCCB
               IDIMAIA = NOCCA*NVIRTA
               IDIMAIB = NOCCB*NVIRTB
               NDET(IRREP)   = IDIMAIA
               NDET(IRREP+8) = IDIMAIB
               NDIMRPA = NDIMRPA + IDIMAIA + IDIMAIB
            END DO
         ELSE
            CALL GETREC(1,'JOBARC','NMPROTON',IONE,NPROT)
            ICHARGE = IFLAGS(h_IFLAGS_charge)
            NELEC   = NPROT-ICHARGE
            NOCCB   = NELEC/2
            NOCCA   = NELEC-NOCCB
            NVIRTA  = NBAS-NOCCA
            NVIRTB  = NBAS-NOCCB
            NDIMRPA = NOCCA*NVIRTA + NOCCB*NVIRTB
         END IF
      END IF
C
      CALL SYMSIZ(NBAS,IRPSZA,IRPSZB,IPKSIZ,MAXMEM,NAOBK,MXLSSZ1,
     &            MXLSSZ2,IRPPMEM,INOUT,NDIMRPA,MAXCOR,IUHF)
C
      MAXMEM=MAXMEM+MOD(MAXMEM,2)
C
      IMEMGET=MAXMEM+IRPPMEM+IFLAGS(44)*IFLAGS(37)
      IMEMGET=IMEMGET+MOD(IMEMGET,2)
cwc1
c more memory is needed if incremental matrices are used and/or 
c Schwarz estimators are used...
      if (dirscf) then
         iTmp = 0
         if (schwrz) then
            call shlask(nsh)
            nsh2 = (nsh*(nsh+1))/2
            iTmp = iTmp + 2*nsh2*IINTFP
         end if
         if (nprocs.gt.1) then
c           buffer for MPI_ALLREDUCE
            iTmp = iTmp + ITRILN(NIRREP+1)*IINTFP
         end if
         if (fdiff) then
            iTmp = iTmp + (1+iuhf)*2*(ITRILN(NIRREP+1)*IINTFP)
         end if
         imemget = imemget + iTmp
      end if
cwc0
C
c OEP SCF block, Igor Schweigert, Sep 2003,
c $Id: vscfoep.FPP,v 1.13 2008/06/06 18:10:43 taube Exp $
C
c An OEP calculation needs more (much more in case of OEP2) memory
c then a regular SCF. I could have modified SYMSIZ to calculate the
c amount of memory needed for different OEP methods, but I want to
c keep the number of modified VSCF files at minimum. So let's ignore
c SYMSIZ's efforts and simply ask for as much memory as specified
c via the MEMORY ASV.

      if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4) then
         if (imemget .lt. iflags (36)) then
            call oep_mem_init(MAXMEM, IRPPMEM, IMEMGET)
         endif
      endif
C
      CALL ACES_MALLOC(IMEMGET,ICORE,I0)
      IF (ICORE(1).EQ.0) THEN
         WRITE(*,*) "@VSCF: Unable to allocate ",IMEMGET,
     &              " words of memory for AO integral sorting."
         CALL ERREX
      END IF
cwc1
c i011 = exchange integrals
c i012 = density matrix compressed to shells
c i013 = MPI_ALLREDUCE buffer
c i021 = alpha density from previous iteration
c i022 =    alpha Fock from previous iteration
c i023 =  beta density from previous iteration
c i024 =     beta Fock from previous iteration
      i011 = i0
      i012 = i0
      i013 = i0
      i021 = i0
      i022 = i0
      i023 = i0
      i024 = i0
      if (dirscf) then
         if (schwrz) then
            i012 = i011 + nsh2*IINTFP
            i0   = i012 + nsh2*IINTFP
         end if
         if (nprocs.gt.1) then
            i013 = i0
            i0   = i013 + ITRILN(NIRREP+1)*IINTFP
         end if
         if (fdiff) then
            i021 = i0
            i022 = i021 + ITRILN(NIRREP+1)*IINTFP
            i0   = i022 + ITRILN(NIRREP+1)*IINTFP
            if (iuhf.eq.1) then
               i023 = i0
               i024 = i023 + ITRILN(NIRREP+1)*IINTFP
               i0   = i024 + ITRILN(NIRREP+1)*IINTFP
            end if
         end if
         imemget = imemget - iTmp
      end if
cwc0

c   o initialize the I/O subsystem for PK chaining
      call aces_io_init(icore,i0,imemget,.true.)
      call aces_init_chemsys
#ifndef _IMPORT_CORR_ORBS
      call aces_io_reset
#endif
c(redundant)      if (INOUT.eq.1) call aces_io_remove(50,'MOINTS')

      MAXCOR=IMEMGET
      I000=I0
c      I010=I000+NBAS+MOD(NBAS,2)
c      I020=I010+NBAS+MOD(NBAS,2)
c      I030=I020+NBAS+MOD(NBAS,2)
c      I040=I030+8*NBAS
c      I050=I040+8*NBAS*IINTFP
c      I060=I050+NBAS*NBAS*IINTFP
c      I070=I060+NBAS+MOD(NBAS,2)
c      I080=I070+NBAS+MOD(NBAS,2)
c      I090=I080+NBAS+MOD(NBAS,2)
c      I100=I090+NBAS+MOD(NBAS,2)
c      IF (I100-I000.GT.MAXMEM) THEN
c        CALL NOMEM('Get SO --> AO transformation matrix','{GETAOSO}'//
c     &             ' <-- VSCF',I100-I000,MAXMEM)
c      END IF
c      CALL GETAOSO(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
c     &             ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
c     &             ICORE(I080),ICORE(I090),NBAS,YESNO2)
C
      CALL TIMER(1)
      IF (IFLAGS(1).GE.10) WRITE(LUOUT,'(a,f5.0,a,/)')
     &   "@VSCF: Initialization and symmetry analysis required ",
     &   TIMENEW," seconds."
C
C DEAVG argument was added on 08/2001 to fix a bug in OOFCOR-DMPSCF
C See comments in DMPSCF.
C 
      IF (INOUT.EQ.2.OR.INOUT.EQ.3) THEN
C NSIZE was unitilized and proper initilization was added, 05/2006,
C Ajitp perera.
C
         CALL GETREC(20,'JOBARC','NAOBASFN',1,NBASX)
         CALL GETREC(20,'JOBARC','NATOMS',1,NATOMS)
         CALL GETREC(-1,'JOBARC','FULLORDR',1,IORDGRP)
         NSIZE = MAX(NBAS,NBASX,NIRREP,NATOMS)
         CALL OOFCOR(ICORE(I000),MAXMEM,IRPPMEM,IBUFLN,NAOBK,MXLSSZ1,
     &               MXLSSZ2,TOL,REPULS,NBAS,IUHF,YESNO2,QC_THRESH,
     &               DAMP,DMPFLG,DMPTOL,DEAVG,SEWARD,IORDGRP,NSIZE)
         GOTO 9889
      END IF
C
      IF (AOFIL) THEN
         I000=I0
         I010=I000
C OPEN IIII INTEGRAL FILE FOR READING 1E INTEGRALS IN GET1E
         if (seward) then
            call gfname('ONE_INT ',fname,ilength)
            open(lunitse,file=fname(1:ilength),form='UNFORMATTED',
     &           access='SEQUENTIAL')
         else
            CALL GFNAME('IIII    ',FNAME,ILENGTH)
            OPEN(LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &           ACCESS='SEQUENTIAL')
         end if
      ELSE
         I000=I0
         I010=I000+IPKSIZ*IINTFP
         I020=I010+IBUFLN*IINTFP
         I030=I020+IBUFLN
         IF (I030-I000.GT.MAXMEM)
     &      cALL NOMEM('Form PK lists',
     &                 '{MKPK2} <-- VSCF',I030-I000,MAXMEM)
         CALL MKPK2(ICORE(I000),ICORE(I010),ICORE(I020),IBUFLN,
     &              IPKSIZ,NBAS,IUHF)
      END IF
C 
C
c OEP SCF block, Igor Schweigert, Sep 2003,
c $Id: vscfoep.FPP,v 1.13 2008/06/06 18:10:43 taube Exp $
c OEP pre-SCF block: <> Check if the OEP method has been invoked,
c set the flags, and notify the user. <> Check fo the known
c incompatibilities. <> Allocate memory to parse the OEP ASV
c list. <> Check if there is enough memory available. <> Check if
c exists and open the ZMAT file. <> Update the OEP control
c variables. <> Close the file. <> Read the number of auxiliary
c basis functions from the OEPINT file. <> In case of an OEP2
c calculation take over the control over convergence, copy the
c threshold into OEPTOL and set the threshold ridiculously small. <>
c Set the iteration type to ' '.
      if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4) then
         call oep_init(petite_list,iuhf,icore(i010),
     &                 (maxmem-i010+i000),loslater)
      endif
C
C PICK UP THE CORE HAMILTONIAN AND OVERLAP MATRICES
C
      I020=I010+ITRILN(NIRREP+1)*IINTFP
      I025=I020+ITRILN(NIRREP+1)*IINTFP+MOD((I020-I000),IINTFP)
      I030=I025+ITRILN(NIRREP+1)*IINTFP
      I040=I030+IBUFLN*IINTFP
      I050=I040+IBUFLN
      I060=I050+NBAS
      I070=I060+ITRILN(NIRREP+1)*IINTFP
      
      IF (I060-I000.GT.MAXMEM)
     &   CALL NOMEM('One electron integrals','{GET1E} <-- VSCF',
     &             I060-I000,MAXMEM)
C
      CALL GET1E(ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &           ICORE(I050),ICORE(I025), ITRILN(NIRREP+1),IBUFLN,
     &           REPULS,NBAS, NIRREP,NBFIRR)
C
C  NOW WE DIAGONLIZE THE OVERLAP MATRIX (STORED STARTING AT I010) AND
C  OBTAIN THE S**(-1/2) ORTHOGONALIZATION MATRIX.  THIS MATRIX IS
C  RETURNED IN PLACE OF THE OVERLAP MATRIX.
C
      I040=I030+MXIRR2*IINTFP
      I050=I040+MXIRR2*IINTFP
      I060=I050+MXIRR2*IINTFP
      I070=I060+MXIRR2*IINTFP
      I080=I070+NBAS*NBAS*IINTFP
      I090=I080+NBAS*IINTFP
      I100=I090+NBAS
      IF (I100-I000.GT.MAXMEM) THEN
         CALL NOMEM('Form S**-1/2','{DOVRLP} <-- VSCF',I100-I000,MAXMEM)
      END IF
C
      CALL DOVRLP(ICORE(I020),ICORE(I025),ICORE(I030),ICORE(I040),
     &            ICORE(I050),ICORE(I060),ICORE(I070),ICORE(I080),
     &            ICORE(I090),ITRILN(NIRREP+1),MXIRR2,NBAS,NLINDEP)
c
c    BEGIN: Find out if where doing a Douglas-Kroll-Hess calc and do it if so.
c
      dkhess  = (iflags2(h_IFLAGS2_dkh_order) .gt. 0)
      if (dkhess) then
         dkhorder=iflags2(h_IFLAGS2_dkh_order)
         if (dkhorder .EQ. 1) DKHZMAX = 1 + 5*NBAS +  8*MXIRR2
         if (dkhorder .EQ. 2) DKHZMAX = 1 + 4*NBAS + 11*MXIRR2
         if (dkhorder .EQ. 3) DKHZMAX = 1 + 4*NBAS + 11*MXIRR2
         if (dkhorder .EQ. 4) DKHZMAX = 1 + 4*NBAS + 21*MXIRR2
         if (dkhorder .EQ. 5) DKHZMAX = 1 + 4*NBAS + 37*MXIRR2

         isize = nbas * (nbas + 1) / 2
         dzero=0.0D0
         call initintgrt(kshf,intnumradpts,dzero)
         call init_ks
         ipnull = setptr(1, F_INTEGER, 1)
         ip010 = setptr(1, F_INTEGER, ibufln)

         iznull = setptr(1, F_REAL, 1)
         iz010 = setptr(1, F_REAL, ibufln)
         iz020 = setptr(1, F_REAL, itriln(nirrep+1))
         iz030 = setptr(1, F_REAL, itriln(nirrep+1))
         iz040 = setptr(1, F_REAL, itriln(nirrep+1))
         iz050 = setptr(1, F_REAL, itriln(nirrep+1))
         iz060 = setptr(1, F_REAL, mxirr2)
         iz070 = setptr(1, F_REAL, mxirr2)
         iz080 = setptr(1, F_REAL, mxirr2)
         iz090 = setptr(1, F_REAL, mxirr2)
         iz100 = setptr(1, F_REAL, isize)
         iz110 = setptr(1, F_REAL, DKHZMAX)
         call dkh__init(dkhorder, nbas, isize,
     +                 ibufln, kscore(ip010),dcore(iz010),
     +                 dcore(iz020), dcore(iz030), dcore(iz040),
     +                 dcore(iz050), dcore(iz060), dcore(iz070),
     +                 dcore(iz080), dcore(iz090), dcore(iz100),
     +                 dcore(iz110), DKHZMAX,
     +
     +                                              icore(i010) )

         call relptr(1, F_REAL, iznull)
         call relptr(1, F_INTEGER, ipnull)
      end if
c
c ...END: Douglas-Kroll-Hess hamiltonian is built, and stored, starting,
c         from icore(i000).
C
C  INITIALIZE THE DENSITY MATRIX TO ZERO.
C
      I040 = I030 + (   (IUHF+1)*ITRILN(NIRREP+1)*IINTFP
     &                + MOD((I030-I000),IINTFP)          )
      CALL ZERO(ICORE(I030),(IUHF+1)*ITRILN(NIRREP+1))
cwc1
c previous-iteration density is initialized to zero, Fock to one-electron
      if (fdiff) then
         call zero (icore(i021),ITRILN(NIRREP+1))
         call zero (icore(i022),ITRILN(NIRREP+1))
         CALL matadd (icore(i010),one,icore(i022),zer,ITRILN(NIRREP+1))
         if (iuhf.ne.0) then
          call zero (icore(i023),ITRILN(NIRREP+1))
          call zero (icore(i024),ITRILN(NIRREP+1))
          CALL matadd (icore(i010),one,icore(i024),zer,ITRILN(NIRREP+1))
         end if
      end if
c calculate exchange integrals
      if (dirscf.and.schwrz) call GAMexc (icore(i011),intg76)
cwc0

      CALL GETREC(20,'JOBARC','NAOBASFN',1,NBASX)
      CALL GETREC(20,'JOBARC','NATOMS',1,NATOMS)
      CALL GETREC(-1,'JOBARC','FULLORDR',1,IORDGRP)
      NSIZE=MAX(NBAS,NBASX,NIRREP,NATOMS)
      I050 = I040 + IINTFP*NBAS*NBAS
      I060 = I050 + IINTFP*MXIRR2
      I070 = I060 + IINTFP*(IUHF+1)*NBAS
      I080 = I070 + IINTFP*MXIRR2
      I090 = I080 + IINTFP*MXIRR2
      I100 = I090 + IINTFP*MXIRR2
      I110 = I100 + IINTFP*NBAS*NBAS
      I120 = I110 + IINTFP*NBAS*NBAS*(IUHF+1)
      I130 = I120 + IINTFP*NBASX*NBAS*(IUHF+1)
      I140 = I130 + IINTFP*NBASX*NBAS
      I150 = I140 + IINTFP*NBASX
      I160 = I150 + IINTFP*NBASX
      I170 = I160 + IINTFP*NBASX
      I180 = I170 + IINTFP*NBAS*(1+IUHF)
      I190 = I180 + IINTFP*NATOMS*3
      I200 = I190 + IINTFP*NATOMS*3
      I210 = I200 + IINTFP*NATOMS
      I220 = I210 + 4*NSIZE*NSIZE        * IINTFP
     &            + MAX(NBASX*10,NATOMS) * IINTFP
     &            + MAX(NBASX+MOD(NBASX,2),IINTFP*3*NATOMS)
     &            + MAX(NBASX+MOD(NBASX,2),NATOMS+MOD(NATOMS,2))
     &            + MAX(IINTFP*NBASX,NATOMS+MOD(NATOMS,2))
     &            + MAX(1080,NATOMS+MOD(NATOMS,2),IINTFP*NBASX*5)
     &            + MAX(2000,NATOMS+MOD(NATOMS,2),NBASX+MOD(NBASX,2),
     &                  NBASX*NBASX*IINTFP)
     &            + MAX(NATOMS*5+MOD(NATOMS,2),IINTFP*NBASX)
     &            + NATOMS*5+MOD(NATOMS,2)
     &            + IORDGRP              * IINTFP
cYAU
c   Messy? We are importing John's modifications for an initial guess and
c rather than check for enough memory when it is needed, all of that stuff
c goes here. Furthermore, ICORE(I210) is still "general" core space, but
c the indexing math must be changed to doubles.
      IF (I220-I000.GT.MAXMEM) THEN
         CALL NOMEM('Initial guess','{INITGES} <-- VSCF',
     &              I220-I000,MAXMEM)
      ELSE
         CALL INITGES(ICORE(I010),ICORE(I020),
     &                ICORE(I030),ICORE(I040),
     &                ICORE(I050),ICORE(I060),
     &                ICORE(I070),ICORE(I080),
     &                ICORE(I090),ICORE(I100),
     &                ICORE(I110),ICORE(I120),
     &                ICORE(I130),ICORE(I140),
     &                ICORE(I150),ICORE(I160),
     &                ICORE(I170),ICORE(I180),
     &                ICORE(I190),ICORE(I200),
     &                ITRILN(NIRREP+1),MXIRR2,
     &                NBAS,NBASX,NATOMS,IUHF,IOLDOCC,
     &                ICORE(I210),(I000+MAXMEM-I220)/IINTFP)
      END IF
C
#ifdef _IMPORT_CORR_ORBS
C This is a highly experiemental work and what this do is
C to read the CC density (or built from natural orbitals) and
C used it subsequently to build a Fock matrix. By choosing
C proper SCF key-words, one can perform various opertaions
C on the resulting Fock matrix (ie. digonalizations). This
c should be kept turned off general purpose ACES II builds.
C Ajith Perera, 05/2009.

C
      CALL GETREC(0, "JOBARC", "EXTORBS", ILENGTH, SCR1)
      CALL GETREC(0, 'JOBARC', 'IP_LIST', ILENGTH, IKEEP)
      CALL GETREC(0, 'JOBARC', 'EA_LIST', ILENGTH, IKEEP)
C
      ILENGTH = 1
      COT  = .TRUE.
      DENS = .TRUE.
      ORBS = .FALSE.
      IF (ILENGTH .GT. 0) THEN
         WRITE(6,*)
         WRITE(6,"(a)") "Correlated Orbitals are being used"
         LDIM1 = ITRILN(NIRREP+1)
         LDIM2 = MXIRR2
         IF (COT .AND. DENS) THEN
            CALL GET_CORR_DENS(ICORE(I030), ICORE(I040), ICORE(I050),
     &                         ICORE(I070), ICORE(I080), ICORE(I090), 
     &                         ICORE(I100), ICORE(I110), ICORE(I120),
     &                         LDIM1, LDIM2, NBAS, NBASX, IUHF)
         ELSE IF (COT .AND. ORBS) THEN
            CALL GETREC(1,'JOBARC','NMPROTON',IONE,NPROT)
            ICHARGE = IFLAGS(h_IFLAGS_charge)
            NELEC   = NPROT-ICHARGE
            NOCCB   = NELEC/2
            NOCCA   = NELEC-NOCCB 
            NVIRTA  = NBAS-NOCCA
C
            I230 = I220 + NOCCA*NOCCA*IINTFP
            I240 = I230 + NVIRTA*NVIRTA*IINTFP
            I250 = I240 + NBAS*IINTFP
            I260 = I250 + NBAS*IINTFP
            I270 = I260 + NOCCA*NOCCA*IINTFP
            I280 = I270 + NOCCA*NOCCA*IINTFP
            I290 = I280 + NVIRTA*NVIRTA*IINTFP
            I300 = I290 + NVIRTA*NVIRTA*IINTFP

            MEMLEFT = (I000 + MAXMEM - I300)/IINTFP
            IF (MEMLEFT .LE. 0) CALL INSMEM("VSCF", I300, MEMLEFT)
C
            CALL GET_CORR_ORBS(ICORE(I030), ICORE(I040), ICORE(I100),
     &                         ICORE(I110), ICORE(I120), ICORE(I130),
     &                         ICORE(I210), ICORE(I070), ICORE(I080),
     &                         ICORE(I090), ICORE(I220), ICORE(I230),
     &                         ICORE(I240), ICORE(I250), ICORE(I260),
     &                         ICORE(I270), ICORE(I280), ICORE(I290),
     &                         ICORE(I300),
     &                         MEMLEFT, LDIM1, LDIM2, NOCCA, NVIRTA,
     &                         NBAS, NBASX, IUHF)
         ENDIF
CSSS         CALL ACES_FIN
CSSS         STOP
C
#ifdef _DEBUG_IMPORT_CORR_ORBS
          CALL DBG_CORR_DENS(ICORE(I030), ICORE(I100), ICORE(I110),
     &                       ICORE(I120), LDIM1, NBAS, NBASX, IUHF)
#endif
       ENDIF 
#endif
c-------------------------------------------------------
      IF (IFLAGS2(h_IFLAGS2_scf_type) .eq.1)  scfks = .TRUE.
      IF (IFLAGS2(h_IFLAGS2_scf_type) .eq.2)  hfdft = .TRUE.
      IF (IFLAGS2(h_IFLAGS2_scf_type) .ge.3)  oepks = .TRUE.
      IF (IFLAGS2(h_IFLAGS2_grad_calc).eq.1)  agrad = .TRUE.
      IF (IFLAGS(h_IFLAGS_props).eq.1)        props = .TRUE.

c prakash for oep and KS mixed 
           if((loslater.or.scfks.or.hfdft))then
           if (.not. dkhess) call initintgrt(kshf,intnumradpts,
     &                                       islaterON)
c lets not turn on the hfdftgrad when the props=first_order. This 
c is added to do HF-DFT like calcs using CC densities. The HFDFT 
c optimizations sould not be combined with property calculations.
c 07/02 Ajith Perera,

           hfdftgrad=hfdft .and. agrad .and. .not. props
           if (hfdftgrad) then
           call putrec(20,'JOBARC','HFDFTGRA',1,1)
           else
            call putrec(20,'JOBARC','HFDFTGRA',1,0)
           end if
           if (.not. dkhess) call init_ks
           call putrec(1,'JOBARC','KSPRINT',iintfp,scfks)
           call setupgrd(max_angpts)
           call getrec(1,'JOBARC','COMPNORB',1,icompnorb)
           ncount = icompnorb
           ic1symmet = 0
           if(CMPGRP .eq.'C1  ') ic1symmet=1
           if (c1symmet .eq. 1) ncount=natoms

           ipnull = setptr(1,F_INTEGER,1)
           iznull = setptr(1, F_REAL,1)
           if(scfks .or.  hfdft) then
           izvalao = setptr(1,F_REAL,nbasx*ncount*intnumradpts
     &                                   *max_angpts)
           izvalgradao = setptr(1,F_REAL,nbasx*ncount*
     &                          intnumradpts*max_angpts*3)

           iztotwt = setptr(1,F_REAL,ncount*intnumradpts 
     &                      *max_angpts)    

           irrtritot= itriln(nirrep+1)*iintfp
           izkspot  = setptr(1, F_REAL,nbasx*nbasx*(iuhf+1))
           iz1      = setptr(1, F_INTEGER,          nbasx*2)
           izscrexc = setptr(1, F_REAL,       nbasx*nbasx*3)
           izscr    = setptr(1,F_REAL,nbasx*nbasx*(iuhf+1))
           izksa    = setptr(1, F_REAL,           irrtritot)
           izksb    = setptr(1, F_REAL,           irrtritot)

           call setnumint(dcore(izvalao),dcore(izvalgradao)
     &                   ,dcore(iztotwt),ncount)

#ifdef _DEBUG
        Print*, "Entering sb_com_mol"
        Print*, "@VSCF.F Iznull, izkspot, iz1, izscrexc, 
     &          izksa, izksb"
        Print*, Iznull, izkspot, iz1, izscrexc, izksa, izksb
#endif
C
c scfexcat in not working has not being incorporated  

           if (scfksexact) then
              call getrec(1,'JOBARC','NAUXBAS',1,numauxbas)
              open(unit=3,file='NUMAUX',status='OLD',form=
     &             'FORMATTED')
            read(3,13) numauxbas
            close(3,status='KEEP')
   13       format(i5)
C
            izscr2   = setptr(1, F_REAL, nbas*nbas)     
            iz2      = setptr(1, F_REAL,nbasx*nbasx*numauxbas)
            izresp   = setptr(1,F_REAL,3*numauxbas*numauxbas)
          end if 
           else
           izvalao = znull                                     
           izvalgradao = znull                         
           iztotwt = znull                                 
           irrtritot= znull                       
           izkspot  = znull                                 
           iz1      = znull                                 
           izscrexc = znull                                 
           izscr    = znull                                
           izksa    = znull                                 
           izksb    = znull                                   
         end if 
           iscre=setptr(1, F_REAL, nbas*nbas)
           iscrq=setptr(1, F_REAL, nbasx*nbas)
           iscrv=setptr(1, F_REAL, nbas*(iuhf+1))
           iscrc=setptr(1,F_REAL, nbas*nbasx*2)
           iscra=setptr(1,F_REAL, nbasx*nbasx)
           iscrocc=setptr(1,F_INTEGER, nbas*2) 
           izend =setptr(1,F_REAL,1)
           ipend =setptr(1,F_INTEGER,1)
      end if
C
      scfkslastiter= .false. 
      scfksiter    = scfks
      If (.not. loslater) then
         slater       = .false.
      else
         slater       = .true.
      endif
C
C This record is read in vdint during HF-DFT gradient calcs. If the
C record is present then vdint will skip the relaxation term (the
C one the involve overlap) since it is computed fully in vksdint.
C Ajith Perera, 02/2008
C
      if (kshf) call putrec(20, "JOBARC", "EXCHANGE", 1, zer)

C-------------------------------------------------------
 
         
      IF (MAXIT.NE.0) THEN
         WRITE(LUOUT,'(T3,68("-"))')
         WRITE(LUOUT,'(T3,"Iteration",9X,"Total Energy",12X,
     &         "Largest Error in FDS-SDF")')
         WRITE(LUOUT,'(T3,68("-"))')
      END IF
C
C  NOW WE DO THE HARD PART:  CONSTRUCT THE FOCK MATRIX.
C
C    I000 = PK LISTS (IN CASE OF AOFIL=.FALSE.)
C    I010 = CORE HAMILTONIAN
C    I020 = ORTHOGONALIZATION MATRIX
C    I030 = OLD DENSITY MATRIX OR MATRICES (UPDATED AT END OF SCF CYCLE)
C    I040 = FOCK MATRIX OR MATRICES
C    I050 = SCRATCH SPACE FOR PK LISTS



C
C  We also iterate over this section.  This is an SCF cycle.
      NOCONV=.TRUE.
      DMAX=0.D0
      ITER=0
      if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4)
     &    call oeptimer ('VSCF: PRE-SCF', 'OFF')
C
 2000 CONTINUE
       ITER=ITER+1
c
c OEP SCF block, Igor Schweigert, Sep 2003,
c $Id: vscfoep.FPP,v 1.13 2008/06/06 18:10:43 taube Exp $
C
       transfer = .false.
       if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4) then
           i070 = i060 + nbas*nbas*iintfp
           i080 = i070 + (iuhf+1)*isqrln(nirrep+1)*iintfp
           if (maxmem-i070+i000 .le. 0) call nomem("do oep iters.",
     &         " oep_dowork <-- vscf", i070-i000, maxmem)
           CALL PUTREC(20,'JOBARC','OCCUPYA0',NIRREP,NOCC(1))
           call oep_dowork(iter, iuhf, nbas, nocc, icore(i070), 
     &                     maxmem-i070+i000, i070-i000,
     &                     icore (i010), icore (i050), 
     &                     icore (i060), icore (i030), icore (i040), 
     &                     etot, MAX_ERROR, itertype, transfer,
     &                     slater,loslater,iSlaterON)
#ifdef _BUILD_KSPOT
        if(transfer) then 
         Call buildKSpot(dcore(izkspot),kscore(iz1),
     &                   dcore(izscrexc),dcore(izscr),
     &                   dcore(izvalao),dcore(izvalgradao),
     &                   dcore(iztotwt),
     &                   intnumradpts,max_angpts,natoms,
     &                   dcore(izksa),ICORE(I040),
     &                   nbasx,0,ITRILN(NIRREP+1))
         end if
#endif
       endif
       if (transfer) go to 5899
C
c scf ks
C
      if (scfksiter .and. num_acc_ks .and. iter .eq. 2) then
         iradtem = int_numradpts
         ipottem = int_defpotgrid
      end if

c num_acc_ks is the switch for the numerical accelerator in the
c KS calculations.
c Also, num_acc_ks prevents the program from
c oscillating around the point where
c the grids are changed. In other words, the SCF routine goes only
c once through the coarse grid.

      if (.not. num_acc_ks) go to 51

      if (scfksiter .and. iter .ge. 2) then
         if (dmax .gt. 0.1d0) then
            int_numradpts  = 30
            int_defpotgrid = 1
        else
             num_acc_ks=.false.
          int_numradpts  = iradtem
          int_defpotgrid = ipottem
         end if
      end if

 51     continue
C
CJDW 4/18/95
C     ----- Piotr's solvent effect stuff. Send problems to him ! -----

C If RFSCF, get the initial H matrix. Audited and fixed silly
C mistakes. I assume that this is correct now eventhough I can't
C exactly reproduce Qchem results. I attribute that to various
C parameters that are being used rather than to a bug. I have
C not tried to reduce mem requirmments (lot of memory being 
C wasted). This repeated call to GETH is not necessary. I 
C wonder why Piotr did not realize the fact that one AO Hamiltonian
C is not suppose to change during iterations.
C Ajith Perera, 06/2008.

      IF (IFLAGS2(118).NE.0) THEN
         I060=I050+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
         I070=I060+IBUFLN*IINTFP
         ITOP=I070+IBUFLN*IINTFP
         CALL GETH(ICORE(I010),ICORE(I050),ICORE(I060),
     &             ITRILN(NIRREP+1),IBUFLN)
c
C adding the reaction field contribution
C
c I060-dipole x
c I070-dipole y
c I080-dipole z
c I090-buf(600)
c I100-ibuf(600)
c I110-scr for expanded density matrix
c I120-scrATCH
C
          I060 = I050
          I070 = I060+NBAS*NBAS*IINTFP
          I080 = I070+NBAS*NBAS*IINTFP
          I090 = I080+NBAS*NBAS*IINTFP
          I100 = I090+ITRILN(NIRREP+1)*IINTFP
          I110 = I100+ITRILN(NIRREP+1)*IINTFP
          I120 = I110+ITRILN(NIRREP+1)*IINTFP
          I130 = I120+ILNBUF*IINTFP
          I140 = I130+ILNBUF*IINTFP
          I150 = I140+NBAS*NBAS*IINTFP
          I160 = I150+MXIRR2*IINTFP
C
          IF (I160-I000.GT.MAXMEM)
     &    CALL NOMEM('Frm rfield',
     &               '{RFIELD2} <-- VSCF',I160-I000,MAXMEM)
          CALL RFIELD2(ICORE(I010),ICORE(I030),ICORE(I060),
     &                 ICORE(I070),ICORE(I080),ICORE(I090),
     &                 ICORE(I100),ICORE(I110),ICORE(I120),
     &                 ICORE(I130),ICORE(I140),ICORE(I150),REPULS,
     &                 ITRILN(NIRREP+1),NBAS,ILNBUF,IUHF,
     &                 MXIRR2,ITER)
      ENDIF
C     ----- End of solvent effect stuff -----

      IF (AOFIL) THEN

C DIRECT CONSTRUCTION OF THE FOCK MATRICES FROM THE AO INTEGRAL FILE
         IF (IUHF.EQ.0) THEN
c         o RHF
            I055=I050 + MOD((I050-I000),IINTFP)
c ICORE(I055) => BUF
c ICORE(I056) => IBUF
c ICORE(I057) => IMAP
            if (petite_list) then
c            o BUF will be the square-packed density matrix. Since this is
c              used piecemeal, we could build each block on the fly and
c              save memory.
c            o IBUF will be the Work array (cf. MOLCAS/scf/ftwo.f)
c            o IBUFLN will measure the maximum free memory of doubles in ICORE
               I056=I055+iintfp*isqrln(nirrep+1)
               imaxbas=nbfirr(iiamax(nirrep,nbfirr,1))
               llt=(imaxbas*(imaxbas+1)/2)
               ineed=(imaxbas*imaxbas)+llt+1
               I057=I056+iintfp*ineed
               if ((I057-I000).gt.maxmem) then
                  call nomem("Make Fock matrix","{MKRHFF} <-- VSCF",
     &                       I057-I000,maxmem)
               else
                  itmp=(maxmem-(I056-I000))/iintfp-1
               end if
            else
               ITMP=IBUFLN
               I056=I055+IBUFLN*IINTFP
               I057=I056+IBUFLN
               IEND=I057+NBAS*NBAS
               IF (IEND-I000.GT.MAXMEM)
     &            CALL NOMEM("Make Fock matrix","{MKRHFF} <-- VSCF",
     &                  IEND-I000,MAXMEM)
            end if
cwc1
c	write (*,*) 'D: current'
c	call dump (icore(i030))
c	write (*,*) 'D: stored'
c	call dump (icore(i021))
        if (fdiff) then
        call matadd (icore(i021),-one,icore(i030),one,ITRILN(NIRREP+1))
        end if
c        write (*,*) 'D: difference'
c        call dump (icore(i030))
        if (dirscf.and.schwrz) then
c prepare the IMAP array in icore(i057); this is also done in MKRHFF,
c but SHLDEN needs this array (BTW, can't it be done once in whole vscf?)
           call mkimap (nirrep,nbas,nbfirr,icore(i057))
c find maxima of shell-pair densities in icore(i012)
           call xSHLDEN
     >       (iuhf,icore(i030),icore(i030),icore(i012),icore(i057),nbas)
        end if
cwc
C-------------------------------------------------------------------
C prakash
C nbasx non symm adapted AO
C iter
C the first iteration must be a HF iteration
C Kohn-sham requires the eigenvector in order to work
C one shoud do i hf iteration
C 0 for iuhf
C ks moudle are such that it need iuhf so I am keeping
C later if i find it is not important or I have means to
C ignore that I will do so
C scfks to provide  the information that it is KS calc
C dcore(iVxc) store Vxc matrix

            CALL MKRHFF(ICORE(I040),
     &                  ICORE(I030),ICORE(I010),
     &                  ICORE(I055),ICORE(I056),
     &                  ITRILN(NIRREP+1),NBAS,NBFIRR,
     &                  ICORE(I057),ITMP,LUINT,.TRUE.,
     &                  nbasx,0,scfksiter,scfksexact,
     &                  scfkslastiter,iter,
     &                  dcore(izkspot),kscore(iz1),
     &                  dcore(izksa),dcore(izscrexc),
     &                  dcore(izscr),dcore(izscr2),dcore(izscreval),
     &                  dcore(izvalao),dcore(izvalgradao),
     &                  dcore(iztotwt),max_angpts,natoms,
     &                  intnumradpts,ncount,kshf)
C small letter means new stuff
C end prakash
C-------------------------------------------------------------------
cwc1
c        write (*,*) 'F: built'
c        call dump (icore(i040))
c        write (*,*) 'F: stored'
c        call dump (icore(i022))
        if (fdiff) then
         call matadd (icore(i021), one,icore(i030),one,ITRILN(NIRREP+1))
c         write (*,*) 'D: summed'
c         call dump (icore(i030))
         call matadd (icore(i022), one,icore(i040),one,ITRILN(NIRREP+1))
         CALL MATADD (icore(i010),-one,icore(i040),one,ITRILN(NIRREP+1))
c         write (*,*) 'F: summed'
c         call dump (icore(i040))
         call matadd (icore(i030),one, icore(i021),zer,ITRILN(NIRREP+1))
         call matadd (icore(i040),one, icore(i022),zer,ITRILN(NIRREP+1))
        end if
cwc0
         ELSE
C         o UHF (AS WELL AS ROHF HIGH SPIN OPEN SHELL)
C
C The following block is added to fix the long unattended
C UHF/FOCK=AO problem (note that the deleterious effect of
C the bug only becomes visible when SYMMETRY=OFF, which explains
C its long survival). 
C Ajith Perera 11/2001
C
            I031=I030+ITRILN(NIRREP+1)*IINTFP
            I032=I031+ITRILN(NIRREP+1)*IINTFP
            I041=I040+ITRILN(NIRREP+1)*IINTFP
            I042=I041+ITRILN(NIRREP+1)*IINTFP
            I050=I042
            I051=I050+MOD((I050-I000),IINTFP)
            I055=I051+ITRILN(NIRREP+1)*IINTFP
c ICORE(I051) => DT
c ICORE(I055) => BUF
c ICORE(I056) => IBUF
c ICORE(I057) => IMAP
            if (petite_list) then
c            o BUF will be the square-packed density matrices.
c            o IBUF will be the Work array (cf. MOLCAS/scf/ftwo.f)
c            o IBUFLN will measure the maximum free memory of doubles in ICORE
               I056=I055+iintfp*isqrln(nirrep+1)*2
               imaxbas=nbfirr(iiamax(nirrep,nbfirr,1))
               llt=(imaxbas*(imaxbas+1)/2)
               ineed=(imaxbas*imaxbas)+llt+1
               I057=I056+iintfp*ineed
               if ((I057-I000).gt.maxmem) then
                  call nomem("Make Fock matrix","{MKUHFF} <-- VSCF",
     &                       I057-I000,maxmem)
               else
                  itmp=(maxmem-(I056-I000))/iintfp-1
               end if
            else
               ITMP=IBUFLN
               I056=I055+IBUFLN*IINTFP
               I057=I056+IBUFLN
               IEND=I057+NBAS*NBAS
               IF (IEND-I000.GT.MAXMEM)
     &            CALL NOMEM("Make Fock matrix","{MKUHFF} <-- VSCF",
     &                  IEND-I000,MAXMEM)
            end if
cwc1
        if (fdiff) then
        call matadd (icore(i021),-one,icore(i030),one,ITRILN(NIRREP+1))
        call matadd (icore(i023),-one,icore(i031),one,ITRILN(NIRREP+1))
        end if
        if (dirscf.and.schwrz) then
c prepare the IMAP array in icore(i057); this is also done in MKUHFF,
c but SHLDEN needs this array (BTW, can't it be done once in whole vscf?)
           call mkimap (nirrep,nbas,nbfirr,icore(i057))
c find maxima of shell-pair densities in icore(i012)
           call xSHLDEN
     >       (iuhf,icore(i030),icore(i031),icore(i012),icore(i057),nbas)
        end if
cwc0
            
            call zero(icore(i040),itriln(nirrep+1))
            CALL MKUHFF(ICORE(I040),ICORE(I041),ICORE(I030),
     &                  ICORE(I031),ICORE(I051),ICORE(I010),
     &                  ICORE(I055),ICORE(I056),
     &                  ITRILN(NIRREP+1),NBAS,NBFIRR,
     &                  ICORE(I057),ITMP,LUINT,.TRUE.,
     &                  nbasx,iuhf,scfksiter,scfksexact,
     &                  scfkslastiter,
     &                  dcore(izkspot),kscore(iz1),
     &                  dcore(izksa),dcore(izksb),dcore(izscrexc),
     &                  dcore(izscr),dcore(izscr2),
     &                  dcore(izvalao),dcore(izvalgradao),
     &                  dcore(iztotwt),max_angpts,natoms,
     &                  intnumradpts,ncount,kshf)

cwc1
        if (fdiff) then
         call matadd (icore(i021), one,icore(i030),one,ITRILN(NIRREP+1))
         call matadd (icore(i022), one,icore(i040),one,ITRILN(NIRREP+1))
         CALL MATADD (icore(i010),-one,icore(i040),one,ITRILN(NIRREP+1))
         call matadd (icore(i030),one, icore(i021),zer,ITRILN(NIRREP+1))
         call matadd (icore(i040),one, icore(i022),zer,ITRILN(NIRREP+1))
         call matadd (icore(i023), one,icore(i031),one,ITRILN(NIRREP+1))
         call matadd (icore(i024), one,icore(i041),one,ITRILN(NIRREP+1))
         CALL MATADD (icore(i010),-one,icore(i041),one,ITRILN(NIRREP+1))
         call matadd (icore(i031),one, icore(i023),zer,ITRILN(NIRREP+1))
         call matadd (icore(i041),one, icore(i024),zer,ITRILN(NIRREP+1))
        end if
cwc0
         END IF
C 
         If (iflags2(h_IFLAGS2_newvrt) .EQ. 1) Then
#ifdef _DEBUG_LVL0
         if (Iter .Eq. 1) 
     &   Write(6,"(a)") "The H + (N-1)/N G is being built"
#endif 

         If (Iuhf .Eq. 0) Then
            Call Make_vnm1_pot(ICORE(I010),ICORE(I040),ICORE(I040),
     &                        Iuhf,ITRILN(NIRREP+1),NOcc)
         Else
            Call Make_vnm1_pot(ICORE(I010),ICORE(I040),ICORE(I041),
     &                        Iuhf, ITRILN(NIRREP+1),NOcc)
         Endif 

         Endif 
C
C Store the Fock matrices in JOBARC. This step is
C required to SCF_MAXCYC=0 option to work correctly.
C 08/07, Ajith Perera.
C
         If (iflags(16) .eq. 0) then

            IADD_SCR3 = I055
            IADD_SCR2 = IADD_SCR3 + NBAS*NBAS*(iuhf+1)*iintfp
            IFBLT_END = IADD_ScR2 + MXIRR2*IINTFP
C
            do ispin=1,(iuhf+1)
               call zero(icore(IADD_SCR3+(ispin-1)*nbas*nbas*iintfp),
     &                    nbas*nbas)
               call mkfull(icore(IADD_SCR3 + (ispin-1)*nbas*nbas*
     &                     iintfp), icore(I040 + (ispin-1)*itriln
     &                     (nirrep + 1)*IINTFP),
     &                     icore(IADD_SCR2), itriln(nirrep + 1),
     &                     mxirr2, nbas, 1)
               label = 'FOCKA   '
               if (ispin .eq. 2) label = 'FOCKB   '
               call putrec(20,'JOBARC', label, nbas*nbas*iintfp,
     &                     icore(iadd_scr3 + (ispin-1)*nbas*nbas*
     &                     iintfp))
            enddo
#ifdef _DEBUG_LVLM
      write(*,*) "The checksum of the new-fock matrices "
      call checksum("FOCK-A", icore(i055), nbas*nbas)
      call checksum("FOCK-B", icore(i055+nbas*nbas*iintfp), nbas*nbas)
CSSS      call output(icore(i055),1,nbas,1,nbas,nbas,nbas,1)
CSSS      call checksum("FOCK-B", icore(i055+nbas*nbas), nbas*nbas)
#endif 
         endif

      ELSE
         I050=I040+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
         I060=I050+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
         I070=I060+MXIRR2*IINTFP
         i080=i070+nbas*nbas*iintfp
         IF (I080-I000.GT.MAXMEM)
     &      CALL NOMEM('Make Fock matrix','{MKFOCK} <-- VSCF',
     &            I080-I000,MAXMEM)
         CALL MKFOCK(ICORE(I000),ICORE(I010),ICORE(I040),ICORE(I030),
     &               ICORE(I050),ICORE(I060),icore(i070),IPKSIZ,
     &               ITRILN(NIRREP+1),MXIRR2,nbas,IUHF)
      END IF
C------transferred from the oep_dowork block-----
 5899 continue
C------

      if (iflags(11).ne.2) then
C Check convergence (RHF/UHF)
        max_error=5.0d0
        i090=i080+nbas*nbas*iintfp
        i100=i090+ISQRLN(NIRREP+1)*IINTFP
        i110=i100+MXIRR2*IINTFP
        i120=i110+mxirr2*iintfp 
        IF (I120-I000.GT.MAXMEM)
     &     CALL NOMEM('Check Conv','{CONVCRIT} <-- VSCF',I110-I000,
     &                MAXMEM)
C----
C i060 location is used in do_oepwork to store eigenvectors. As
C a consequence this block needed to be replaced (along with the
C arguments to concrit). Ajith Perera, 07/2009. 
CSSS        i060=i050+nbas*nbas*iintfp
CSSS        I070=I060+ISQRLN(NIRREP+1)*IINTFP
CSSS        I080=I070+MXIRR2*IINTFP
CSSS        i090=i080+mxirr2*iintfp 
CSSS        IF (I090-I000.GT.MAXMEM)
CSSS     &     CALL NOMEM('Check Conv','{CONVCRIT} <-- VSCF',I090-I000,
CSSS     &                MAXMEM)
C---
        call convcrit(icore(i040),icore(i030),icore(i080),icore(i090),
     &                icore(i100),icore(i110),itriln(nirrep+1),nbas,
     &                iuhf,max_error)
      else if (iter.eq.1) then
        max_error=5.0d0
      endif
C
C  Now compute the total energy.
c prakash oep
c since we calculated energy in oep routine thus we skip mkenerg
c transfer true meaning OEP
c end prakash
C
      if (.not. transfer) then

      I060=I050+ITRILN(NIRREP+1)*IINTFP
      I070=I060+MXIRR2*IINTFP
      I080=I070+MXIRR2*IINTFP
      IF (I080-I000.GT.MAXMEM)
     &   CALL NOMEM('Form energy','{MKENER} <-- VSCF',I080-I000,MAXMEM)  
      CALL MKENER(ICORE(I010),ICORE(I030),ICORE(I040),ICORE(I050),
     &            ICORE(I060),ICORE(I070),ITRILN(NIRREP+1),MXIRR2,
     &            REPULS,mAX_ERROR,ITER-1,IUHF,ETOT,0,NOCONV,scfksiter,
     &            scfkslastiter,EDIFF)
      end if
c end transfer

      if (max_error.lt.tol) then
        if (qc_phase) then
          qc_phase=.false.
          iflags(h_IFLAGS_scf_extrapo)=0
        else
          noconv=.false.
        endif
       else if (max_error.lt. 1.0E-04 .AND. EDIFF .LT. 10E-12 .AND.
     &  NLINDEP .NE. 0)then
          Write(*,*)
          Write(*,"(A,A)") " Warning: The density did not fall",
     &          " below the desired thershold and the density"
          Write(*,"(A,A)") " is converged only to 10E-05. However, the",
     &                     " the energy is converged to 10E-12."
          Write(*,"(A,A)") " The calculation is allowed to proceed", 
     &                     " since linear dependencis are found."
          Write(*,"(A,A)") " caution is urged for using refrence",
     &                     " for post-HF methods."
          noconv=.false.
       endif 
c
      if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4)
     &    call oeptimer ('SCFITER: REST', 'ON')
c
C  Do Davidson or static damping.
cTODO: Should we add .AND..NOT.QC_PHASE ?
      IF (DMPFLG)
     &    CALL DMPSCF(ICORE(I040),ICORE(I050),ITRILN(NIRREP+1),ETOT,
     &                ITER,IUHF,DAMP,DMPFLG,DMPTOL,DEAVG)
     
      IF (ROHFMO   .AND.IFLAGS(11).EQ.2) THEN
         I041=I040+ITRILN(NIRREP+1)*IINTFP
         I031=I030+ITRILN(NIRREP+1)*IINTFP
         I051=I050+IINTFP*ITRILN(NIRREP+1)
         I052=I051+IINTFP*ITRILN(NIRREP+1)
         I053=I052+IINTFP*ISQRLN(NIRREP+1)
         I054=I053+IINTFP*ISQRLN(NIRREP+1)
         IEND=I054+IINTFP*ISQRLN(NIRREP+1)
         CALL MKROHFF(ICORE(I040),ICORE(I041),
     &                ICORE(I050),ICORE(I051),
     &                ICORE(I052),ICORE(I053),
     &                ICORE(I054),NBFIRR,NOCC,ITRILN(NIRREP+1),
     &                ISQRLN(NIRREP+1))
        call convcrit(icore(i040),icore(i030),icore(iend),icore(i052),
     &                icore(i053),icore(i054),itriln(nirrep+1),nbas,
     &                iuhf,max_error)
        if (max_error.lt.tol) noconv=.false.
      END IF
C  IF THIS CALCULATION USES DIIS, THEN GO AHEAD AND DO IT.
      IF (IFLAGS(10).NE.0.AND..NOT.QC_PHASE.AND.ITER.GT.1) THEN
         NDIMA=MAX(MXIRR2,(IORDER+1)*(IORDER+2))
C
C         NDIMB=MAX(MXIRR2,IORDER+1)
C To check the stability of the RPP matrix.
C
         NDIMB=MAX(MXIRR2,(IORDER+1)*(IORDER+1))
         I060=I050+NBAS*NBAS*IINTFP
         I070=I060+(IUHF+1)*ISQRLN(NIRREP+1)*IINTFP
         I080=I070+MXIRR2*IINTFP
         I090=I080+MXIRR2*IINTFP
         I100=I090+MXIRR2*IINTFP
         I110=I100+MAX(MXIRR2,(IORDER+1)*(IORDER+2))*IINTFP
c        I120=I110+MAX(MXIRR2,IORDER+1)*IINTFP
c To match NDIMB, above
         I120=I110+MAX(MXIRR2,(IORDER+1)*(IORDER+1))*IINTFP
         I130=I120+MXIRR2*IINTFP
         IF (I130-I000.GT.MAXMEM)
     &      CALL NOMEM('Convergence acceleration','{DORPP} <-- VSCF',
     &            I130-I000,MAXMEM)
         CALL DORPP(ICORE(I000+MAXMEM),ICORE(I040),ICORE(I030),
     &              ICORE(I050),ICORE(I060),ICORE(I070),ICORE(I080),
     &              ICORE(I090),ICORE(I100),ICORE(I110),ICORE(I120),
     &              ITRILN(NIRREP+1),ISQRLN(NIRREP+1),IRPPMEM,
     &              MXIRR2,NDIMA,NDIMB,NBAS,ITER-1,IUHF,TOL,NOCONV)
      ENDIF

C  Now we have fun.  Let's do the rest of the SCF iteration now.
C  Pass all this wonderful stuff into the next routine.
C
C     I000 = PJ lists (only in case AOFIL=.FALSE.)
C     I010 = Core Hamiltonian
C     I020 = Orthogonalization matrix
C     I030 = Density matrix
C     I040 = Fock matrix
C     I050 = Vector for orbital energies (both alpha and beta)
C     I060 = Orbital eigenvectors
C     I070 = Scratch density matrix
C     I080 = Scratch matrix
C     I090 = Scratch matrix
C     I100 = Scratch matrix
C     I110 = Scratch matrix

      I060=I050+(IUHF+1)*NBAS*IINTFP
      I070=I060+(IUHF+1)*ISQRLN(NIRREP+1)*IINTFP
      I080=I070+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      I090=I080+MAX(2*NBAS,NBAS*NBAS,MXIRR2)*IINTFP
      I100=I090+MAX(2*NBAS,NBAS*NBAS,MXIRR2)*IINTFP
      I110=I100+MAX(2*NBAS,MXIRR2)*IINTFP
      I120=I110+MAX(2*NBAS,MXIRR2)*IINTFP

C Two additional scratch arrays needed if orbital locking has been
C requested (C(OLD)T * S * C is monitored). If it has not been
C requested, we just use ICORE(I100) and ICORE(I110) for these arrays.

      IF (GSSLOK) THEN
         I130 = I120 + NBAS * NBAS * IINTFP
         I140 = I130 + NBAS * NBAS * IINTFP
      ELSE
         I120 = I100
         I130 = I110
         I140 = I130 + MAX(2*NBAS,MXIRR2)*IINTFP
      END IF

      IF (I140-I000.GT.MAXMEM)


     &   CALL NOMEM('SCF iteration','{SCFIT} <-- VSCF',I140-I000,MAXMEM)

C ADD IN MEMORY FOR QCSCF
C FIRST COMPUTE # OF ORBITAL ROTATION COEFFICIENTS (SINGLE EXCITATIONS)
C IN EACH IRREP BY SYMMETRY AND SPIN SO WE KNOW HOW MUCH SPACE WE NEED.
C TOTAL # IS NDIM AND STORE BY IRREP AND SPIN IN NDET ARRAY
      NDIM=0
      DO J=1,NIRREP
         NOCCA=NOCC(J)
         NOCCB=NOCC(J+8)
         NVIRTA=NBFIRR(J)-NOCCA
         NVIRTB=NBFIRR(J)-NOCCB
         NDET(J)=NOCCA*NVIRTA
         NDET(J+8)=NOCCB*NVIRTB
         NDIM=NDIM+NDET(J)+NDET(J+8)
      END DO

      I150=I140+IINTFP*NDIM ! FVEC
      I160=I150+IINTFP*NDIM ! XTRIAL
      I170=I160+IINTFP*NBAS*NBAS ! CVECS
      I180=I170+IINTFP*NBAS*NBAS ! SOVER
      I190=I180+IINTFP*NDIM ! RESIDUE
      I200=I190+IINTFP*NDIM ! D1
      I210=I200+IINTFP*NDIM ! D2
      I220=I210+IINTFP*NDIM ! D3
      I230=I220+IINTFP*NDIM ! D4
      I240=I230+IINTFP*ISQRLN(NIRREP+1)*(IUHF+1) ! RAA
      I250=I240+IINTFP*ISQRLN(NIRREP+1)*(IUHF+1) ! G1AA
      I260=I250+IINTFP*ISQRLN(NIRREP+1)*(IUHF+1) ! G2AA

      IF (I260-I000.GT.MAXMEM.AND.QC_PHASE)
     &   CALL NOMEM('QC acceleration','{QCSCF} <-- SCFIT <-- VSCF',
     &              I260-I000,MAXMEM)
          
C The primary responsiblity of checking convergence has been transfered 
C to RPP when RPP is on. The SCFIT still handle the cases where RPP is
C turned off by the user or quadratic convergence procedure is used. 
C Ajith Perera 12/07

      IF (IFLAGS(16).NE.0) THEN
         CALL SCFIT(ICORE(I010),ICORE(I020),ICORE(I030),ICORE(I040),
     &              ICORE(I050),ICORE(I060),ICORE(I070),ICORE(I080),
     &              ICORE(I090),ICORE(I100),ICORE(I110),ICORE(I120),
     &              ICORE(I130),REPULS,     ITRILN(NIRREP+1),MXIRR2,
     &              ISQRLN(NIRREP+1),NOCONV,TOL,        ITER,
     &              MAX_ERROR,            IUHF,       NBAS,
     &              ICORE(I140),ICORE(I150),ICORE(I160),ICORE(I170),
     &              ICORE(I180),ICORE(I190),ICORE(I200),ICORE(I210),
     &              ICORE(I220),ICORE(I230),ICORE(I240),ICORE(I250),
     &              NPRVOCC,    IPREVOCC,
     &              NDIM,       NDET,       QC_THRESH,  QC_PHASE,
     &          dcore(iscrq),scfks.or.hfdft.or.loslater,dcore(iscre),
     &              dcore(iscrc),dcore(iscrv),dcore(iscra),nbasx,
     &              kscore(iscrocc))
      ELSE
         CALL RELOAD(ICORE(I060),ICORE(I050),ICORE(I080),ICORE(I090),
     &               ICORE(I100),ISQRLN(NIRREP+1),NBAS,IUHF)
         noconv=.false.
      END IF
C
      if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4) then
         call oeptimer ('SCFITER: REST', 'OFF')
         call oeptimer ('VSCF: SCFITER', 'OFF')
      endif
c
       if (scfks) then
          call relptr(1,F_REAL,izend)
          call relptr(1,F_INTEGER,ipend)
       end if
C
      IF (NOCONV.AND.ITER.LT.MAXIT) GOTO 2000
C
      if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4) 
     &    call oeptimer ('VSCF: POSTSCF', 'ON')
c
cjp modified diagnostics
      IF (NOCONV.or.iflags(16).eq.0) THEN
cjp scf process skipped when maxit==0
         if (iflags(16).eq.0) then
            write(luout,'(/,a,/)') "@VSCF: SCF has been omitted."
         else
            if (iflags(16).eq.1) then
               write(luout,'(/,a,/)') "@VSCF: one iteration performed."
            else
               WRITE(LUOUT,'(/,a,/)') "@VSCF: SCF has not converged!!!"
               KILLFLG(1)=1
            end if
         end if
      ELSE
         WRITE(LUOUT,'(/,a,/)') "@VSCF: SCF has converged."
      END IF
C
c OEP SCF block, Igor Schweigert, Sep 2003,
c $Id: vscfoep.FPP,v 1.13 2008/06/06 18:10:43 taube Exp $
C
      if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4) then
         call oep_post(iuhf, nbas, ibufln,luint,nocc, icore(i070),
     &                 maxmem-i070+i000, icore (i050), 
     &                 icore (i030), etot, dmax, itertype)
         call putrec(20, 'JOBARC', 'SCFEVCA0', nbas*nbas*iintfp, 
     &               Icore(i060))
         CALL MKRHFF(ICORE(I040),
     &               ICORE(I030),ICORE(I010),
     &               ICORE(I055),ICORE(I056),
     &               ITRILN(NIRREP+1),NBAS,NBFIRR,ICORE(I057),
     &               ITMP,LUINT,.TRUE.,
     &               nbasx,0,scfksiter,scfksexact,scfkslastiter,
     &               iter,
     &               dcore(izkspot),kscore(iz1),
     &               dcore(izksa),dcore(izscrexc),
     &               dcore(izscr),dcore(izscr2),dcore(izscreval),
     &               dcore(izvalao),dcore(izvalgradao),
     &               dcore(iztotwt),max_angpts,natoms,
     &               intnumradpts,ncount,kshf)
         call getrec(20, 'JOBARC', 'SCFEVLA0', nbas*iintfp, 
     &               icore(i050))
         call getrec(20, 'JOBARC', 'SCFEVCA0', nbas*nbas*iintfp, 
     &               Icore(i060))
         goto 5999
      endif
C
      if (.not. noconv .AND. scfks)then
          scfksiter=scfks
          scfkslastiter=.true.
      else if(.not. noconv .AND. hfdft) then
         scfksiter=hfdftgrad
         scfkslastiter=.false.
      end if
C
C  Construct the final fock matrix and the final energy.
      I080=I070+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
      I090=I080+MXIRR2*IINTFP
      i100=i090+nbas*nbas*iintfp
      IF (I100-I000.GT.MAXMEM)
     &   CALL NOMEM('Final Fock matrix','{MKFOCK} <-- VSCF',I100-I000,
     &         MAXMEM)
      IF (AOFIL) THEN
         IF (IUHF.EQ.0) THEN
            I055=I070 + MOD((I070-I000),IINTFP)
c ICORE(I055) => BUF
c ICORE(I056) => IBUF
c ICORE(I057) => IMAP
            if (petite_list) then
               I056=I055+iintfp*isqrln(nirrep+1)
               llt=(imaxbas*(imaxbas+1)/2)
               ineed=(imaxbas*imaxbas)+llt+1
               I057=I056+iintfp*ineed
               if ((I057-I000).gt.maxmem) then
                  call nomem("Make Fock matrix","{MKRHFF} <-- VSCF",
     &                       I057-I000,maxmem)
               else
                  itmp=(maxmem-(I056-I000))/iintfp-1
               end if
            else
               ITMP=IBUFLN
               I056=I055+IBUFLN*IINTFP
               I057=I056+IBUFLN
               IEND=I057+NBAS*NBAS
               IF (IEND-I000.GT.MAXMEM)
     &            CALL NOMEM("Make Fock matrix","{MKRHFF} <-- VSCF",
     &                  IEND-I000,MAXMEM)
            end if
cwc1
        if (dirscf.and.schwrz) then
CSSS           call mkimap (nirrep,nbas,nbfirr,icore(i057))
CSSS           call xSHLDEN
CSSS     >       (iuhf,icore(i030),icore(i030),icore(i012),icore(i057),nbas)
        end if
cwc0
            CALL MKRHFF(ICORE(I040),
     &                  ICORE(I030),ICORE(I010),
     &                  ICORE(I055),ICORE(I056),
     &                  ITRILN(NIRREP+1),NBAS,NBFIRR,ICORE(I057),
     &                  ITMP,LUINT,.TRUE.,
     &                  nbasx,0,scfksiter,scfksexact,scfkslastiter,
     &                  iter,
     &                  dcore(izkspot),kscore(iz1),
     &                  dcore(izksa),dcore(izscrexc),
     &                  dcore(izscr),dcore(izscr2),dcore(izscreval),
     &                  dcore(izvalao),dcore(izvalgradao),
     &                  dcore(iztotwt),max_angpts,natoms,
     &                  intnumradpts,ncount,kshf)
C
         ELSE IF (IUHF.EQ.1) THEN
C         o UHF
C The following block is added to fix the long unattended
C UHF/FOCK=AO problem (note that the deleterious effect of
C the bug only becomes visible when SYMMETRY=OFF, which explains  
C its long survival).
C Ajith Perera 11/2001
C
            I031=I030+ITRILN(NIRREP+1)*IINTFP
            I032=I031+ITRILN(NIRREP+1)*IINTFP
            I041=I040+ITRILN(NIRREP+1)*IINTFP
            I042=I041+ITRILN(NIRREP+1)*IINTFP
C
CAjith Perera 11/2001
C 
            I051=I070+MOD((I050-I000),IINTFP)
            I055=I051+ITRILN(NIRREP+1)*IINTFP
c ICORE(I051) => DT
c ICORE(I055) => BUF
c ICORE(I056) => IBUF
c ICORE(I057) => IMAP
            if (petite_list) then
c            o BUF will be the square-packed density matrices.
c            o IBUF will be the Work array (cf. MOLCAS/scf/ftwo.f)
c            o IBUFLN will measure the maximum free memory of doubles in ICORE
               I056=I055+iintfp*isqrln(nirrep+1)*2
               llt=(imaxbas*(imaxbas+1)/2)
               ineed=(imaxbas*imaxbas)+llt+1
               I057=I056+iintfp*ineed
               if ((I057-I000).gt.maxmem) then
                  call nomem("Make Fock matrix","{MKUHFF} <-- VSCF",
     &                       I057-I000,maxmem)
               else
                  itmp=(maxmem-(I056-I000))/iintfp-1
               end if
            else
               ITMP=IBUFLN
               I056=I055+IBUFLN*IINTFP
               I057=I056+IBUFLN
               IEND=I057+NBAS*NBAS
               IF (IEND-I000.GT.MAXMEM) CALL
     &            NOMEM("Make Fock matrix","{MKUHFF} <-- VSCF",
     &                  IEND-I000,MAXMEM)
            end if
cwc1
        if (dirscf.and.schwrz) then
           call mkimap (nirrep,nbas,nbfirr,icore(i057))
           call xSHLDEN
     >       (iuhf,icore(i030),icore(i031),icore(i012),icore(i057),nbas)
        end if
cwc0
            CALL MKUHFF(ICORE(I040),ICORE(I041),ICORE(I030),
     &                  ICORE(I031),ICORE(I051),ICORE(I010),
     &                  ICORE(I055),ICORE(I056),
     &                  ITRILN(NIRREP+1),NBAS,NBFIRR,
     &                  ICORE(I057),IBUFLN,LUINT,.TRUE.,
     &                  nbasx,iuhf,scfksiter,scfksexact,
     &                  scfkslastiter,
     &                  dcore(izkspot),kscore(iz1),
     &                  dcore(izksa),dcore(izksb),dcore(izscrexc),
     &                  dcore(izscr),dcore(izscr2),
     &                  dcore(izvalao),dcore(izvalgradao),
     &                  dcore(iztotwt),max_angpts,natoms,
     &                  intnumradpts,ncount,kshf)
C
         END IF
C
C For HF-DFT gradients modified Fock matrix. Prakash Verma and Ajith Perera
C 042008
C
         if (hfdftgrad) then
C
            do is=1,iuhf+1
              do ir=1,nirrep
                if(nbfirr(ir).ne.0)then
                ioccp(is)=nocc((is-1)*8+ir)
                ivrtp(is)=nbfirr(ir)-ioccp(is)
                ioccvrt(is)=ioccvrt(is)+ioccp(is)*ivrtp(is)
                end if
              end do
            end do

            istart_scr  = i055
            istart_fock = istart_scr  + nbas*nbas*iintfp
            iscrtmp     = istart_fock + mxirr2*iintfp
            ict         = iscrtmp     + mxirr2*iintfp
            iftia       = ict         + mxirr2*iintfp
            iftib       = iftia       +ioccvrt(1)*iintfp
            ifblt_end   = iftib  +ioccvrt(iuhf+1)*iintfp

            call zero(icore(istart_scr), nbas*nbas*iintfp)
            call mkfull(icore(istart_scr), icore(I040), 
     &                  icore(istart_fock), itriln(nirrep + 1),
     &                  mxirr2, nbas, 1)

            call putrec(1,'JOBARC','HFKSVDAA',nbas*nbas*
     &                  iintfp,icore(istart_scr))
C
            if(iuhf .eq.1) then
              call zero(icore(istart_scr), nbas*nbas*iintfp)
              call mkfull(icore(istart_scr), icore(I041),
     &                  icore(istart_fock), itriln(nirrep + 1),
     &                  mxirr2, nbas, 1)
            end if
C
            call putrec(1,'JOBARC','HFKSVDBB',nbas*nbas*iintfp
     &                  ,icore(istart_scr))

            call blt_hfksFAI(icore(i060),iuhf,isqrln(nirrep+1),
     &                      icore(i040),icore(istart_fock),
     &                      mxirr2,itriln(nirrep+1),
     &                      icore(iscrtmp),icore(ict),icore(iftia),
     &                      icore(iftib),ioccvrt(1),ioccvrt(2),
     &                      nocc)
         endif
C 
C Store the Fock matrices in JOBARC. This step is
C required to SCF_MAXCYC=0 option to work correctly.
C 08/07, Ajith Perera.
C
         If (iflags(16) .eq. 0) then
            IADD_SCR3 = I055
            IADD_SCR2 = IADD_SCR3 + NBAS*NBAS*(iuhf+1)*iintfp
            IFBLT_END = IADD_ScR2 + MXIRR2*IINTFP
C
            do ispin=1,(iuhf+1)
               call zero(icore(IADD_SCR3+(ispin-1)*nbas*nbas*iintfp),
     &                    nbas*nbas)
               call mkfull(icore(IADD_SCR3 + (ispin-1)*nbas*nbas*
     &                     iintfp), icore(I040 + (ispin-1)*itriln
     &                     (nirrep + 1)*IINTFP),
     &                     icore(IADD_SCR2), itriln(nirrep + 1),
     &                     mxirr2, nbas, 1)
               label = 'FOCKA   '
               if (ispin .eq. 2) label = 'FOCKB   '
               call putrec(20,'JOBARC', label, nbas*nbas*iintfp,
     &                     icore(iadd_scr3 + (ispin-1)*nbas*nbas*
     &                     iintfp))
            enddo
C
C If the QRHF flags are set, adjust the new template occupations
C with the QRHF orbital options. Otherwise, we do it here. 
C Note that the OCCNUM file should correspond to the closed shell
C (just like standard QRHF). 09/2015, Ajith perera. 

C The functions implemented in this procedure can be alternatively acomplished 
C by using the OCCNUM file and no longer in use. 07/2021. 
CSSS            If (IFLAGS(77).EQ.0) THEN
CSSS               CALL GETREC(-1,'JOBARC','OSCALC  ',1,IOS)
CSSS               CALL TEMPLATE_ORBS(ICORE(I010),ICORE(I030),ICORE(I040),
CSSS     &                           ICORE(I050),ICORE(I060),
CSSS     &                           ITRILN(NIRREP+1),
CSSS     &                           ISQRLN(NIRREP+1),NBAS,IUHF,IOS)
CSSS            Endif 
#ifdef _DEBUG_LVLM
      write(*,*) "The checksum of the new-fock matrices "
      call checksum("FOCK-A", icore(i055), nbas*nbas)
      call checksum("FOCK-B", icore(i055+nbas*nbas*iintfp), nbas*nbas)
CSSS      call output(icore(i055),1,nbas,1,nbas,nbas,nbas,1)
CSSS      call checksum("FOCK-B", icore(i055+nbas*nbas), nbas*nbas)
#endif
         endif

      ELSE
         CALL MKFOCK(ICORE(I000),ICORE(I010),ICORE(I040),ICORE(I030),
     &               ICORE(I070),ICORE(I080),icore(i090),IPKSIZ,
     &               ITRILN(NIRREP+1),MXIRR2,nbas,IUHF)
      END IF

c   o Now that we are done calling MKxHFF, we need to close the MOLCAS
c     2-electron integral file.
      if (petite_list) then
         itmp=-1
         call clsord(itmp,0)
         if (itmp.ne.0) then
            write(*,*) "@VSCF: ERROR - A problem was encountered ",
     &                 "while attempting to close"
            write(*,*) "       the MOLCAS 2-electron integral file."
            call errex
         end if
      end if

C GET ORBITAL ENERGIES FOR ROHF
      IF (IFLAGS(11).EQ.2) THEN
         I041=I040+ITRILN(NIRREP+1)*IINTFP
         I051=I070
         I052=I051+ISQRLN(NIRREP+1)*IINTFP
         I053=I052+ISQRLN(NIRREP+1)*IINTFP
         I054=I053+ISQRLN(NIRREP+1)*IINTFP
         CALL GETEVAL(ICORE(I040),ICORE(I041),ICORE(I050),
     &                ICORE(I051),ICORE(I052),ICORE(I053),
     &                ITRILN(NIRREP+1),ISQRLN(NIRREP+1),NBFIRR)
      END IF

      if (iflags(11).ne.2) then
        i080=i070+nbas*nbas*iintfp
        I090=I080+ISQRLN(NIRREP+1)*IINTFP
        I100=I090+MXIRR2*IINTFP
        i110=i100+mxirr2*iintfp 
        IF (I110-I000.GT.MAXMEM)
     &     CALL NOMEM('Check Conv','{CONVCRIT} <-- VSCF',I100-I000,
     &              MAXMEM)
        max_error=5.0d0
        call convcrit(icore(i040),icore(i030),icore(i070),icore(i090),
     &                icore(i100),icore(i110),itriln(nirrep+1),nbas,
     &                iuhf,max_error)

      endif

      I080=I070+ITRILN(NIRREP+1)*IINTFP
      I090=I080+MXIRR2*IINTFP
      I100=I090+MXIRR2*IINTFP
      IF (I100-I000.GT.MAXMEM) CALL
     &   NOMEM('Final energy','{MKENER} <-- VSCF',I100-I000,MAXMEM)
      CALL MKENER(ICORE(I010),ICORE(I030),ICORE(I040),ICORE(I070),
     &            ICORE(I080),ICORE(I090),ITRILN(NIRREP+1),MXIRR2,
     &            REPULS,MAX_ERROR,ITER,IUHF,ETOT,0,NOCONV,scfksiter,
     &            scfkslastiter,EDIFF)

       if(ikk.eq.2) goto 3008
C
C------- transferred from oep_post block----
 5999 continue
C---------
C
C  Print out the final occupations, if changed.
      CHGOCC=.FALSE.
      DO I=1,16
         IF (NOCC(I).NE.IOLDOCC(I)) CHGOCC=.TRUE.
      END DO
      IF (CHGOCC) THEN
         WRITE(LUOUT,'(/,a,/)') "@VSCF: Final occupancies:"
         WRITE(LUOUT,'(T8,"   Alpha population by irrep: ",8(I3,2X))')
     &      (NOCC(I),I=1,NIRREP)
         WRITE(LUOUT,'(T8,"    Beta population by irrep: ",8(I3,2X))')
     &      (NOCC(8+I),I=1,NIRREP)
         WRITE(LUOUT,*)
      END IF
 3008 continue
C  For open-shell singlet calculations checks phase of orbitals
C
C
C     Definition of CHKPHASE is provisionally changed. It is assumed
C     that IOS is set in ALL open-shell singlet calculations, which
C     means that REF=TWODET should have been set. Previously OSPHASE
C     was called in ALL UHF calculations with equal alpha and beta
C     occupations.
C
C     Incidentally, there seems to be a potential problem here for
C     open-shell singlet calculations when FOCK=AO, since .TRUE. implies
C     that PK formation is done in core.
C
c& new lines
      IDELOCC = 0
      DO I=1, NIRREP
         IDELOCC = IDELOCC + IABS(NOCC(I) - NOCC(8+I))
      END DO
      CALL GETREC(-1,'JOBARC','OSCALC  ',1,IOS)
cjdw  CHKPHASE = (IOS .NE. 0) .OR. ((IUHF .GE. 1) .AND.
cjdw &   (IDELOCC .EQ. 0))
      CHKPHASE = (IOS.NE.0)
      IF (CHKPHASE) THEN
     
         I080=I070+NBAS*NBAS*IINTFP
         I090=I080+MXIRR2*IINTFP
         I100=I090+MXIRR2*IINTFP
         MXMEM=MAXMEM-(I100-I000)
         IF (I100-I000.GT.MAXMEM) CALL
     &      NOMEM('Rephase orbitals','{OSPHASE} <-- VSCF',I100-I000,
     &            MAXMEM)

CSSS         CALL OSPHASE(ICORE(I060),ICORE(I070),ICORE(I080),
CSSS     &                ICORE(I090),ISQRLN(NIRREP+1),MXIRR2,
CSSS     &                ITRILN(NIRREP+1),NBAS,ICORE(I100),MXMEM,IUHF,
CSSS     &                ICORE(I030),ICORE(I000),ICORE(I010),ICORE(I040),
CSSS     &                IPKSIZ,.FALSE.)

            CALL OSPHASE(ICORE(I000),ICORE(I010),ICORE(I030),
     &              ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I100),
     &              MAXMEM-(IEND-I000),ISQRLN(NIRREP+1),MXIRR2,
     &              NBAS,IPKSIZ,REPULS,ETOT,IUHF,LUINT,
     &              IOS,ICORE(I070),ICORE(I080),ICORE(I090),.FALSE.,
     &              nbasx,scfksiter,scfksexact,scfkslastiter,
     &              dcore(izkspot),kscore(iz1),
     &              dcore(izksa),dcore(izksb),dcore(izscrexc),
     &              dcore(izscr),dcore(izscr2),
     &              dcore(izvalao),dcore(izvalgradao),
     &              dcore(iztotwt),max_angpts,natoms,
     &              intnumradpts,ncount,kshf)
C
      END IF
c& end

C  Now print out the eigenvectors, if the user wants them.
C
CJDW 11/27/96. Stop SO eigenvector printing for the time being. We can
C              print them in more useful AO basis with EVCAO2.
C
c      IF (SCFPRNT.EQ.0) THEN
c        WRITE(LUOUT,5210)
c 5210   FORMAT(T3,'@VSCF-I, Eigenvector printing suppressed.',/)
c      ELSE
c        I080=I070+NBAS+MOD(NBAS,2)
c        I090=I080+NBAS+MOD(NBAS,2)
c        IF (I090-I000.GT.MAXMEM) THEN
c          CALL NOMEM('Dump eigenvectors','{OUTEVC} <-- VSCF',I090-I000,
c     &               MAXMEM)
c        END IF
cC
c        CALL OUTEVC(ICORE(I050),ICORE(I060),
c     &              NBAS,ISQRLN(NIRREP+1),IUHF,AOFIL)
c      END IF
C
C  Now calculate the value of S**2, if UHF.
      IF (IUHF.EQ.1) THEN
         I080=I070+NBAS*NBAS*IINTFP
         I090=I080+MXIRR2*IINTFP
         I100=I090+MXIRR2*IINTFP
         IF (I100-I000.GT.MAXMEM) CALL
     &      NOMEM('Calculate S**2','{S2CALC} <-- VSCF',I100-I000,
     &            MAXMEM)
         CALL S2CALC(ICORE(I060),ICORE(I070),ICORE(I080),
     &               ICORE(I090),ISQRLN(NIRREP+1),MXIRR2,NBAS)
      END IF

C Dump the SCF occupations now for subsequent steps of an optimization.
      CALL PUTREC(20,'JOBARC','OCCUPYA0',NIRREP,NOCC(1))
      IF (IUHF.NE.0)
     &   CALL PUTREC(20,'JOBARC','OCCUPYB0',NIRREP,NOCC(9))

CMN
CMN  transform occ and virt orbitals to v(N-1) potential if
CMN  requested
CMN
      IF (iflags2(h_IFLAGS2_newvrt) .EQ. 2) THEN
         I080=I070+MAX(2*NIRREP,NBAS+MOD(NBAS,2))
         I090=I080+5*ISQRLN(NIRREP+1)*IINTFP
         I100=I090+ITRILN(NIRREP+1)*IINTFP
         IF (i100-I000.GT.MAXMEM) CALL
     &      NOMEM('NEW VIRTUALS','{DRNEWVRT} <-- VSCF',i100-I000,
     &            MAXMEM)
         CALL DRNEWVRT(icore(i010),icore(i020),ICORE(I030),ICORE(I040),
     &                 ICORE(I050),
     &                 ICORE(I060),ICORE(I070),ICORE(I080),ICORE(I090),
     &                 ITRILN(NIRREP+1),
     &                 ISQRLN(NIRREP+1),MXIRR2,NBAS,IUHF)

C  For open-shell singlet calculations checks phase of orbitals
c& new lines
         IF (CHKPHASE) THEN
            I080=I070+NBAS*NBAS*IINTFP
            I090=I080+MXIRR2*IINTFP
            I100=I090+MXIRR2*IINTFP
            MXMEM=MAXMEM-(I100-I000)
            IF (I100-I000.GT.MAXMEM) CALL
     &         NOMEM('Rephase orbitals','{OSPHASE} <-- VSCF',I100-I000,
     &               MAXMEM)

CSSS            CALL OSPHASE(ICORE(I060),ICORE(I070),ICORE(I080),
CSSS     &                   ICORE(I090),ISQRLN(NIRREP+1),MXIRR2,
CSSS     &                   ITRILN(NIRREP+1),NBAS,ICORE(I100),MXMEM,IUHF,
CSSS     &                   ICORE(I030),ICORE(I000),ICORE(I010),
CSSS     &                   ICORE(I040),IPKSIZ,.TRUE.)

            CALL OSPHASE(ICORE(I000),ICORE(I010),ICORE(I030),
     &              ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I100),
     &              MAXMEM-(IEND-I000),ISQRLN(NIRREP+1),MXIRR2,
     &              NBAS,IPKSIZ,REPULS,ETOT,IUHF,
     &              IOS,ICORE(I070),ICORE(I080),ICORE(I090),.TRUE.,
     &              nbasx,scfksiter,scfksexact,scfkslastiter,
     &              dcore(izkspot),kscore(iz1),
     &              dcore(izksa),dcore(izksb),dcore(izscrexc),
     &              dcore(izscr),dcore(izscr2),
     &              dcore(izvalao),dcore(izvalgradao),
     &              dcore(iztotwt),max_angpts,natoms,
     &              intnumradpts,ncount,kshf)

         END IF
c& end

         IF (FIRSTSCF .AND. (SCFPRNT.GE.1)) THEN
            CALL OUTEVC(ICORE(I050),ICORE(I060),
     &                  NBAS,ISQRLN(NIRREP+1),IUHF,AOFIL)
         END IF
      END IF
CMN END
C  If this is a QRHF calculation, set it up. If th

      IF (IFLAGS(77).NE.0) THEN
        CALL QRHFIT(ICORE(I000),ICORE(I010),ICORE(I030),ICORE(I040),
     &          ICORE(I050),ICORE(I060),ICORE(I070),
     &          MAXMEM-(I060-I000),ITRILN(NIRREP+1),
     &          ISQRLN(NIRREP+1),NBAS,IPKSIZ,REPULS,ETOT,IUHF,IOS,
     &          nbasx,scfksiter,scfksexact,scfkslastiter,
     &          dcore(izkspot),kscore(iz1),
     &          dcore(izksa),dcore(izksb),dcore(izscrexc),
     &          dcore(izscr),dcore(izscr2),
     &          dcore(izvalao),dcore(izvalgradao),
     &          dcore(iztotwt),max_angpts,natoms,
     &          intnumradpts,ncount,kshf)
      Endif

C Include section for UNO-REF
C
C     I000 = PJ lists (only in case AOFIL=.FALSE.)
C     I010 = Core Hamiltonian
C     I020 = Orthogonalization matrix
C     I030 = Density matrix
C     I040 = Fock matrix
C     I050 = Vector for orbital energies (both alpha and beta)
C     I060 = Orbital eigenvectors
C     I070 = Scratch density matrix
C
      LDIM2 = ISQRLN(NIRREP+1)
      IDENSNEW=I070 
      IDENSSCR=IDENSNEW+2*LDIM2*IINTFP
      IDENSFRAC=IDENSSCR+2*LDIM2*IINTFP
      ISCRA=IDENSFRAC+2*LDIM2*IINTFP
      ISCRB=ISCRA+2*LDIM2*IINTFP
      IEND=ISCRB+2*LDIM2*IINTFP
      IF (iEND-I000.GT.MAXMEM) CALL
     &   NOMEM('Dump to JOBARC','{DMPJOB} <-- VSCF',iEND-I000,
     &         MAXMEM)
c
cmn write recored nsmican to jobarc. Will be overwritten if needed
c
      if (iflags(16) .ne. 0) then
         i = 0
         call putrec(20,'JOBARC','NSMICAN ', 1, i)
      endif
c
      IF(IFLAGS2(148).NE.0) THEN
C
C Test here changed to be based on makerhf instead of multiplicity.
C Also added is uno_makerhf_light. This do very basic UNO steps
C Ajith Perera, 08/2015.
C 
        if (iflags2(140) .eq. 0) then
        CALL UNO_REF(ICORE(I000),ICORE(I010),ICORE(I030),ICORE(I040),
     &              ICORE(I050),ICORE(I060),ICORE(IEND),
     &              MAXMEM-(IEND-I000),ITRILN(NIRREP+1),
     &              ISQRLN(NIRREP+1),NBAS,IPKSIZ,REPULS,ETOT,IUHF,
     &              LUINT,
     $        ICORE(I020), ICORE(IDENSNEW), ICORE(IDENSSCR),
     $        ICORE(ISCRA), ICORE(ISCRB), ICORE(IDENSFRAC),
     &          nbasx,scfksiter,scfksexact,scfkslastiter,
     &          dcore(izkspot),kscore(iz1),
     &          dcore(izksa),dcore(izksb),dcore(izscrexc),
     &          dcore(izscr),dcore(izscr2),
     &          dcore(izvalao),dcore(izvalgradao),
     &          dcore(iztotwt),max_angpts,natoms,
     &          intnumradpts,ncount,kshf)
        else
         IF (PROGRAM .EQ. "MN_A3") THEN
            CALL UNO_MAKERHF(
     $              ICORE(I000),ICORE(I010),ICORE(I030),ICORE(I040),
     &              ICORE(I050),ICORE(I060),ICORE(IEND),
     &              MAXMEM-(IEND-I000),ITRILN(NIRREP+1),
     &              ISQRLN(NIRREP+1),NBAS,IPKSIZ,REPULS,ETOT,IUHF,
     &              LUINT,
     $              ICORE(I020), ICORE(IDENSNEW), ICORE(IDENSSCR),
     $              ICORE(ISCRA), ICORE(ISCRB), ICORE(IDENSFRAC),
     &              nbasx,scfksiter,scfksexact,scfkslastiter,
     &              dcore(izkspot),kscore(iz1),
     &              dcore(izksa),dcore(izksb),dcore(izscrexc),
     &              dcore(izscr),dcore(izscr2),
     &              dcore(izvalao),dcore(izvalgradao),
     &              dcore(iztotwt),max_angpts,natoms,
     &              intnumradpts,ncount,kshf)
          ELSE
            CALL UNO_MAKERHF_LIGHT(
     $              ICORE(I000),ICORE(I010),ICORE(I030),ICORE(I040),
     &              ICORE(I050),ICORE(I060),ICORE(IEND),
     &              MAXMEM-(IEND-I000),ITRILN(NIRREP+1),
     &              ISQRLN(NIRREP+1),NBAS,IPKSIZ,REPULS,ETOT,IUHF,
     &              LUINT,
     $              ICORE(I020), ICORE(IDENSNEW), 
     $              ICORE(ISCRA), ICORE(ISCRB), ICORE(IDENSSCR),
     &              nbasx,scfksiter,scfksexact,scfkslastiter,
     &              dcore(izkspot),kscore(iz1),
     &              dcore(izksa),dcore(izksb),dcore(izscrexc),
     &              dcore(izscr),dcore(izscr2),
     &              dcore(izvalao),dcore(izvalgradao),
     &              dcore(iztotwt),max_angpts,natoms,
     &              intnumradpts,ncount,kshf)
          ENDIF 

        endif
      ENDIF
C
C  Now dump the stuff to JOBARC.
      CALL GETREC(20,'JOBARC','NAOBASFN',1,NBASX)
      I080=I070+NBAS*NBAS*IINTFP
      I090=I080+MXIRR2*IINTFP
      I100=I090+MAX(2*NIRREP,NBAS+MOD(NBAS,2))
      I110=I100+4*ISQRLN(NIRREP+1)*IINTFP
      I120=I110+NBAS*NBAS*IINTFP
C---Added to create no symmetry AO basis density; 06/2011, Ajith Perera
      I130=I120+NBASX*NBASX*IINTFP
      I140=I130+NBASX*NBAS*IINTFP
      I150=I140+NBASX*NBAX*IINTFP
      IF (i150-I000.GT.MAXMEM) CALL
     &   NOMEM('Dump to JOBARC','{DMPJOB} <-- VSCF',i150-I000,
     &         MAXMEM)
      CALL DMPJOB(icore(i010),icore(i020),ICORE(I030),ICORE(I040),
     &            ICORE(I050),
     &            ICORE(I060),ICORE(I070),ICORE(I080),ICORE(I090),
     &            ICORE(I100),ICORE(I110), 
     &            ICORE(I120), ICORE(I130), ICORE(I140), 
     &            ITRILN(NIRREP+1),
     &            ISQRLN(NIRREP+1),MXIRR2,NBAS,NBASX, ETOT,IUHF,SEWARD)

C Now form the eigenvectors in the original AO basis (not the SO basis).
C EVCAO2 is called regardless of print level. Logic in EVCAO2 decides
C whether vectors are to be printed.
C
CJDW 11/7/96.
C     EVCAO has not been called for ages. The new routine EVCAO2
C     is used to put MOs in ZMAT ordered AO basis (NBASX*NBAS) and print
C     them out. EVCAO2 also writes a formatted file which will soon be
C     used for initial guesses.
c      I080=I070+NBAS*NBAS*IINTFP
c      I090=I080+NBAS*NBAS*IINTFP
c      I100=I090+NBAS*NBAS*IINTFP
c      IF (I100-I000.GT.MAXMEM) THEN
c        CALL NOMEM('Eigenvectors (SO --> AO)','{EVCAO} <-- VSCF',
c     &             I100-I000,MAXMEM)
c      END IF
cc      CALL EVCAO(ICORE(I070),ICORE(I080),ICORE(I090),NBAS,IUHF)
C Moved to before Dmpjob Call, 06/2011, Ajith Perera
CSSS      CALL GETREC(20,'JOBARC','NAOBASFN',1,NBASX)

      CALL GETREC(20,'JOBARC','NATOMS',1,NATOMS)
      I010 = I000 + IINTFP*NBAS *NBAS
      I020 = I010 + IINTFP*NBASX*NBAS
      I030 = I020 + IINTFP*NBASX*NBAS
      I160 = I030 + IINTFP*NBAS *2
      I040 = I160 + IINTFP*3*NATOMS
      I050 = I040 +        NATOMS
      I060 = I050 +        NBASX
      I070 = I060 +        NBASX
      I080 = I070 +        NATOMS
      I090 = I080 +        NATOMS
      I100 = I090 +        NATOMS
      I110 = I100 +        NATOMS
      I140 = I110 +        NATOMS
      I150 = I140 +        NBASX
      I170 = I150 +        NBAS*2
      IF (I170-I000.GT.MAXMEM) THEN
         CALL NOMEM('Eigenvectors (SO --> AO)','{EVCAO2} <-- VSCF',
     &              I170-I000,MAXMEM)
      ELSE
         CALL EVCAO2(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     &               ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
     &               ICORE(I080),ICORE(I090),ICORE(I100),ICORE(I110),
     &               ICORE(I140),NOCC,       ICORE(I150),ICORE(I160),
     &               NBASX,      NBAS,       NATOMS,     IUHF)
      END IF

c  Now we do some fun stuff.  We determine the irreps of each of
c  the eigenvectors.  If this is a finite difference calculation,
c  then don't try to do this; it won't work (calculation in C1 while
c  full irrep stuff in true symmetry of system).
c
c      CALL GETREC(-1,'JOBARC','JODAFLAG',I16,ICORE(I000))
c      IF ((ICORE(I000+8).NE.2) .AND. FIRSTSCF) THEN
      CALL IRRORB_SHELL(ICORE(I000),MAXMEM,NBAS,IUHF)
      CALL PRJDENS(ICORE(I000),MAXMEM,NBAS,IUHF)

C dump the reordered eigenvalues picking them up from JOBARC
      I010 = I000 + (IUHF+1)*NBAS*IINTFP
      I020 = I010 + NBAS+MOD(NBAS,2)
      I030 = I020 + NBAS+MOD(NBAS,2)
      IF (I030-I000.GT.MAXMEM) THEN
         CALL NOMEM('Dump eigenvalues','{EVLOUT} <-- VSCF',
     &              I030-I000,MAXMEM)
      ELSE
         if (firstscf) then
            CALL EVLOUT(ICORE(I000),ICORE(I010),ICORE(I020),
     &                  IREPS,NBAS,NIRREP,IUHF,NBFIRR)
         end if
      END IF

C The SA-AO occupied-occupied block of the Fock matrix is need for the 
C gradients. Built it here. This is need for SCF/Correlated gradients
C with Alaska code. 07/2000, Ajith Perera
C Lets built these records anyway in order to get the density plot
C cords work. Nothing done here should intefere with anything else.
C 08/2015, Ajith Perera.

CSSS      if (seward) then
         postscf = (iflags(2).ne.0)
         I010 = I000 + IINTFP*NBAS*NBAS
         I020 = I010 + IINTFP*NBAS*NBAS
         I030 = I020 + IINTFP*NBAS*NBAS
         I040 = I030 + IINTFP*NBAS*NBAS*2
         I050 = I040 + IINTFP*NBAS*NBAS
         IF ((I050-I000).GT.MAXMEM) THEN
            CALL NOMEM('Eff. Fock','{BLT_EFCTVFOCK} <-- VSCF',
     &                 I050-I000,MAXMEM)
         ELSE
            CALL BLT_EFCTVFOCK(ICORE(I000), ICORE(I010), ICORE(I020),
     &                         ICORE(I030), ICORE(I040), NBAS, NOCC,
     &                         NIRREP, NBFIRR, IUHF, POSTSCF)
         END IF
CSSS      end if
C
C  Now for the really hard stuff . . . the MULLIKEN POPULATION
C  ANALYSIS!!!
C
c     I010=I000+NBAS*NBAS*IINTFP
c     I020=I010+NBAS*NBAS*IINTFP
c     I030=I020+NBAS*NBAS*IINTFP
c     I040=I030+(IUHF+1)*NBAS*IINTFP
c     I050=I040+NBAS+MOD(NBAS,2)
c     I060=I050+NBAS+MOD(NBAS,2)
c     I070=I060+NBAS+MOD(NBAS,2)
c     I080=I070+NBAS+MOD(NBAS,2)
c     I090=I080+NBAS*IINTFP
c     I100=I090+NBAS+MOD(NBAS,2)
c     I110=I100+NBAS+MOD(NBAS,2)
c     I120=I110+NBAS+MOD(NBAS,2)
c     IF (I120-I000.GT.MAXMEM) THEN
c       CALL NOMEM('Mulliken population analysis','{MULPOP} <-- VSCF',
c    &             I120-I000,MAXMEM)
c     END IF
C
c      CALL MULPOP(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
c     &            ICORE(I040),ICORE(I050),ICORE(I060),ICORE(I070),
c     &            ICORE(I080),ICORE(I090),ICORE(I100),ICORE(I110),
c     &            NBAS,IUHF)
c
c for mrcc/aces3 calculations reorder eigenvectors according to
c active labels
c
        if (iflags2(132) .eq. 3 .and. iflags(54) .ne. 0
     $       .and. iflags(16) .ne. 0) then
c
           ieigvec = i000
           ieigval = ieigvec + nbas*nbas*iintfp
           iscr = ieigval+nbas*iintfp
           nscr = 3*nbasx*nbasx + 2 * nbas
           iend = iscr + nscr*iintfp
           if ((iend - i000) .gt. maxmem)
     &      call nomem('mrcc prep','{prep_mrcc} <-- vscf',
     &             iend-i000,maxmem)
c
           call prep_mrcc(nbas, nbasx, icore(ieigvec),
     $          icore(ieigval), icore(iscr), nscr,
     $          nirrep, nocc, nbfirr, iintfp)
c
        endif
C 
c 9889 IF (INOUT.NE.1) CLOSE(MOINTS,STATUS='DELETE')
 9889 IF (INOUT.NE.1) call aces_io_remove(50,'MOINTS')
      inquire(unit=LUINT,opened=bTmp)
      if (bTmp) CLOSE(LUINT,STATUS='KEEP')
C
      DO KILLNUM=1,3
         IF (KILLFLG(KILLNUM).NE.0) THEN
            TERMINATE=.TRUE.
            IF (KILLNUM.EQ.1) THEN
               WRITE(LUOUT,'(a,i3,a,/)')
     &            "@VSCF: SCF failed to converge in ",MAXIT," cycles."
            ELSE IF (KILLNUM.EQ.2) THEN
               WRITE(LUOUT,'(a,/)')
     &            "@VSCF: Linear dependencies in basis."
               WRITE(LUOUT,'(a,/)')
     &            "Only single-point calculations allowed."
            ELSE IF (KILLNUM.EQ.3) THEN
               WRITE(LUOUT,'(a,/)')
     &            "@VSCF: Wavefunction is not symmetric."
            END IF
         END IF
      END DO
C
      if (iflags2 (153).eq.3 .or. iflags2 (153).eq.4) then
         call oeptimer ('VSCF: POSTSCF', 'OFF')
         call oeptimer ('VSCF', 'OFF')
         call oeptimer ('', 'REPORT')
      endif
c
cjp for bwread we will need to know if iuhf BEFORE crapsi is called!
      open(unit=2,file='IUHF',form='unformatted',access='sequential')
      write(2)iuhf
      close(2)
cjp
cwc1
      close(99)
      call MPW_Finalize(iErr)
      if (iErr.ne.MPI_SUCCESS) then
         print *, '@VSCF: Error in MPI_Finalize'
         call aces_exit(1)
      end if
cwc0

      IF (TERMINATE) THEN
C Even though we're about to die, we may want to restart from
C this point, so make sure all the stuff in DMPJOB actually makes it
C to disk.
         call aces_fin
         CALL ERREX
         STOP 1
      ELSE
         CALL ENDSCF
      END IF
      call callstack_term
      END

