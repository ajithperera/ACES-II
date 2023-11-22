
C VECTORIZED CC/MBPT PROGRAM
C
C WRITTEN BY J.F. STANTON, J. GAUSS, J.D. WATTS AND W.J. LAUDERDALE
C TD-CC EXTENSION ADDED BY P.G. SZALAY (ROUTINES WITH 'OS')
cjp MR-BW-CCSD extension added by Jiri Pittner (1998-99)

      PROGRAM VCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISTSZ,RLECYC
      INTEGER POP(8,2),VRT(8,2)
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD,LINCC
      LOGICAL CICALC,NONHF,TRIPIT,TRIPNI,TRIPNI1,T3STOR,PRESNT,INIT,
     &        DORESET,UCC,RESTART,BRUECK,BRKCNV,READT
      LOGICAL ROHF4,ITRFLG,OSFLAG,INITL1
      LOGICAL CIS,EOM
      DIMENSION ECORR(3)
      COMMON / / ICORE(1)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON /TRIPLES/TRIPNI,TRIPNI1,TRIPIT,T3STOR
      COMMON /ENERGY/ ENERGY(500,2),IXTRLE(500)
      COMMON /LINEAR / LINCC,CICALC
      COMMON /NHFREF/ NONHF
      COMMON /ROHF/ ROHF4,ITRFLG
      COMMON /CORENG/ ELAST
      COMMON /BRUECKNER/ BRUECK
      COMMON /SYMLOC/ ISYMOFF(8,8,25)
      COMMON /SYM/    POP,VRT,NT(2),NFMI(2),NFEA(2)
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     &                NRECS,LUAUX
      COMMON /T3IOOF/ IJKPOS(8,8,8,2),IJKLEN(36,8,4),IJKOFF(36,8,4),
     &                NCOMB(4)
      COMMON /EXCITE/ CIS,EOM
cjp
#include "bwcc.com"
#ifdef MPI
#include "mpif.h"
#endif
cjp signal handler need acces to ncycle
      integer ncycle
      common /ncycle/ncycle
cjp
C
      EQUIVALENCE (IFLAGS(2),METHOD)
c
cjp parallelization setup
       call bwparinit
cjp for signal processing
       call installsig
cjp
C
 1000 FORMAT(/,77('-'),/,32X,'Entering xvcc',/,77('-'),/)
      WRITE(6,1000)
C
cjp
      call bwread(.true.,iuhf)
cjp
C
cjp!!! not necessary to repeat, in fact multiple call of popcorn wastes memory
cjp MOINTS will be just opened in crapsi, need such one which will exist
      iref = myproc
cjp
C
      CALL CRAPSI(ICORE,IUHF,0)
C
cjp do not allow crazy attempts
      if (isbwcc.and.
     &    (method.ne.8.and.method.ne.10.and.method.ne.22)) then
         print *,
     &      'MR-BW-CC available only within CCD and CCSD approximations'
         call aces_exit(1)
      end if
cjp we assume that record numbers are the same in all MO files
cjp this could be done parallelized, but is negligible
      if (nref.gt.nproc) then
         do iref = 1, nref
            call storemoio
         end do
      end if
      if (isbwcc.and.(.not.masik)) call bwprep(nocco,nvrto,iuhf)
      if (isbwcc) nonhf=.true.
cjp
C
      CALL ZERO(ENERGY,1000)
      CALL IZERO(IXTRLE,500)
      CALL SETMET
      CALL SETLOG(RESTART)
      ICONVG=1
C
cjp
      if (isbwcc) then
         do i = 1, nref
            ibwconvg(i) = 1
         end do
      end if
cjp
C
      CALL GETREC(-1,'JOBARC','OSCALC',1,IOS)
      OSFLAG = (IOS.NE.0)
      BRUECK = (IFLAGS(22).EQ.1)
      ITRFLG = .TRUE.
      ROHF4  = .FALSE.

c     INITIALIZATION FOR OPENSHELL SINGLET CALCULATION
      IF (OSFLAG) CALL OSINIC(ICORE,ICRSIZ,IUHF,I0)
C
c Nevin added to insure maxcor is aligned from the top as well as the bottom
c      MAXCOR=ICRSIZ
      MAXCOR=ICRSIZ-iand(icrsiz,1)
C
cjp
      if (.not.isbwcc) nref = 1
      if (.not.isbwcc) ibwpass = 0
      do 1234 iref = 1, nref
cjp parallelization
      if (mod(iref-1,nproc).ne.myproc-1) goto 1234
      if (nref.gt.nproc) call reopenmo
C
C INITIALIZE DOUBLES LISTS
      IF (METHOD.GT.1) CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
C
C  INITIALIZE LISTS FOR UCC CALCULATION
      IF (UCC) CALL UCCLST(ICORE(I0),MAXCOR,IUHF)
C
C INITIALIZE SINGLE LISTS
      IF (METHOD.GE.3.OR.NONHF) CALL INSING(METHOD,IUHF,RESTART)
C
C-----------------------------------------------------------------------
C     Create and initialize triples lists if these are needed.
C-----------------------------------------------------------------------
      T3STOR = .FALSE.
      METACT = 0
      IF (TRIPNI.OR.TRIPNI1.OR.TRIPIT) THEN
         T3STOR = METHOD.EQ.12.OR.
     &            METHOD.EQ.17.OR.
     &            METHOD.EQ.18.OR.
     &           (METHOD.GE.26.AND.METHOD.LE.31).OR.
     &           (TRIPIT.AND.IFLAGS(38).GT.0.AND.IFLAGS(39).EQ.0)
c 3/14/97     &   .OR.
c 3/14/97     &        (METHOD.EQ.14.AND.IFLAGS(87).EQ.3).OR.
c 3/14/97     &        (METHOD.EQ.16.AND.IFLAGS(87).EQ.3).OR.
c 3/14/97     &        (METHOD.EQ.33.AND.IFLAGS(87).EQ.3).OR.
c 3/14/97     &        (METHOD.EQ.34.AND.IFLAGS(87).EQ.3)
c11/ 3/96     &        (METHOD.EQ.34.AND.IFLAGS(87).EQ.3).OR.
c11/ 3/96     &        (METHOD.EQ.22.AND.IFLAGS(87).EQ.3 .AND.
c11/ 3/96     &           (IFLAGS2(124).EQ.4 .OR. IFLAGS2(124).EQ. 6 .OR.
c11/ 3/96     &            IFLAGS2(124).EQ.8 .OR. IFLAGS2(124).EQ.10) )
         IF (METHOD.EQ.29.OR.METHOD.EQ.30) THEN
            METACT = METHOD
            METHOD = 18
         END IF
         IF (T3STOR) THEN
 1010       FORMAT('  @VCC-I, Triples lists will be created. ')
            WRITE(6,1010)
C           PRESNT indicates whether file already exists or not.
C           At this stage in a calculation it should not exist.
            PRESNT = .FALSE.
            INIT   = .TRUE.
            IRREPX = 1
            NLIST  = 8
            IF (IFLAGS(87).EQ.3.AND.IFLAGS(91).GT.0) THEN
               IF (METHOD.EQ.13) NLIST = 12
               IF (METHOD.EQ.14) NLIST = 12
               IF (METHOD.EQ.16) NLIST = 12
               IF (METHOD.EQ.33) NLIST = 12
               IF (METHOD.EQ.34) NLIST = 12
               IF (METHOD.EQ.18) NLIST = 16
               IF (METHOD.EQ.22) NLIST =  4
            END IF
            CALL AUXIOI(PRESNT,INIT,IRREPX,NLIST,IUHF)
c        END IF (T3STOR)
         END IF
c     END IF (TRIPNI.OR.TRIPNI1.OR.TRIPIT)
      END IF
C-----------------------------------------------------------------------

C SET UP SOME INFORMATION BEFORE CALCULATION.
C
CJDW  3/20/96 IEVERY is set equal to IFLAGS(13), so that TAMP_SUM
C             keyword works.
C      IEVERY = 5
      IEVERY = IFLAGS(13)
cjp debug...
      if (isbwcc.and.bwgossip) ievery=1
cjp
C
      IF (CC) THEN
         NCYCLE = IFLAGS(7)
         NKEEP  = IFLAGS(12)
         ICONTL = IFLAGS(4)
      END IF
cjp debug...
      IF (isbwcc.or.(IFLAGS(1).GE.10)) THEN
         WRITE(6,801)
  801    FORMAT(T3,' Initial T amplitudes: ')
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
         WRITE(6,800)
  800    FORMAT(T3,' Correlation energies computed from initial ',
     &          'T amplitudes: ')
      END IF
      if (isbwcc) then
         call cmpeng(icore(i0),maxcor,43,0,ecorr,edummy1,edummy2,iuhf,1)
         energy(1,1) = 0.d0
         energy(1,2) = 0.d0
      else
         CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,
     &               ENERGY(1,1),ENERGY(1,2),IUHF,1)
      end if
      ELAST = ENERGY(1,2)
C
C BLOCK OF CODE FOR FINITE ORDER CALCULATIONS.
C
      if (nref.gt.nproc) then
         call aces_cache_flush
         call storemoio
      end if

cjp end iref loop 1
1234  continue

c ----------------------------------------------------------------------

      IF (.NOT.CC) THEN
CMN
CJDW 10/24/95. Block of MN's code.
C
C GENERATE INTERMEDIATES FOR MBPT(2) BASED EOM CALCULATION
         IF (IFLAGS2(117).EQ.2) THEN
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL INSING(10,IUHF,RESTART)
            CALL SETLST(ICORE(I0),MAXCOR,IUHF)
            CALL INITIN(ICORE(I0),MAXCOR,IUHF)
            CALL INITSN(ICORE(I0),MAXCOR,IUHF)
            CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,2)
            CALL INIT2(IUHF)
C
C  PUT CONVERGED AMPLITUDES AS INITIAL GUESS FOR SOLVING THE
C  LAMBDA EQUATIONS ON THE GAMLAM FILE
C
            DO ISPIN = 3, 3-2*IUHF, -1
               CALL MVSYM(ICORE(I0),MAXCOR,IUHF,1,43+ISPIN,143+ISPIN)
            END DO
            initl1 = .false.
            if (initl1) then
               DO ISPIN=1,IUHF+1
                  CALL UPDMOI(1,NT(ISPIN),ISPIN,190,0,0)
                  CALL UPDMOI(1,NT(ISPIN),ISPIN+2,190,0,0)
               END DO
C              ZERO THE L1 LIST (ALL CYCLES)
               CALL ZERLST(ICORE,NT(1),1,1,1,190)
               IF (IUHF.NE.0) CALL ZERLST(ICORE,NT(2),1,1,2,190)
            end if
C
            IF (IUHF.EQ.0) CALL RESET(ICORE(I0),MAXCOR,IUHF)
            IF (IFLAGS(35).NE.0) CALL ACES_AUXCACHE_FLUSH
            WRITE(6,1020)
            call aces_fin
            call bwparend
            call aces_exit(0)
c        END IF (IFLAGS2(117).EQ.2)
         END IF
CMN END

C CALL SPECIAL ROUTINE FOR NONHF PERTURBATION THEORY
         IF (NONHF) THEN
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL ROHFPT(ICORE(I0),MAXCOR,IUHF)
C IN CASE OF MBPT BASED ON BRUECKNER ORBITALS, GET NEW ORBITALS
C BY CALLING BRUECKIT
            IF (BRUECK) THEN
               write(*,5100)
5100           FORMAT('   Residual T1 amplitudes ')
               CALL AMPT1(ICORE(I0),MAXCOR,IUHF)
               CALL BRUECKIT(ICORE(I0),MAXCOR,IUHF,BRKCNV)
            END IF
C FOR GRADIENT CALCULATION, INITIALIZE ``GAMLAM-LISTS'' HERE
            CALL INCOR(I0,IUHF)
            MAXCOR = ICRSIZ
C
            IF (GRAD.AND.METHOD.GT.1) THEN
C
C FOR ROHF-MBPT(3) GRADIENRS, COPY T[2] TO L[2] LISTS
C
C  T[1] LISTS 44,45,46   T[2] LISTS 61,62,63
C  L[1] IS  T[1]         L[2] LISTS 143,144,145
C
C   SINGLES : T[1]: 1,2;90   T[2]:  3,4;90
C             L[1] IS T[1]   L[2]:  1,2;190
C
               CALL INIT2(IUHF)
               CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,.TRUE.)
            END IF
C
C HF CASES
C
         ELSE IF (MBPT3) THEN
            CALL SETLST(ICORE(I0),MAXCOR,IUHF)
            CALL INCOR(I0,IUHF)
            MAXCOR = ICRSIZ
            CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,0)
            CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)
            CALL NEWT2(ICORE(I0),MAXCOR,IUHF)
            CALL CMPENG(ICORE(I0),MAXCOR,60,0,ECORR,ENERGY(2,1),
     &                  ENERGY(2,2), IUHF,0)
            IF (GRAD) CALL FORMT2(ICORE(I0),MAXCOR,IUHF,43,60)
         ELSE IF (MBPT4) THEN
C
C FOURTH-ORDER LOGIC.  COMPUTE D3 AND D4 ON FIRST PASS, THEN
C  GET QUADS.
C
            CALL SETLST(ICORE(I0),MAXCOR,IUHF)
            CALL INCOR(I0,IUHF)
            MAXCOR=ICRSIZ
            CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,0)
            CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)
            CALL NEWT2(ICORE(I0),MAXCOR,IUHF)
            CALL GETD4(ICORE(I0),MAXCOR,ENERGY(3,1),IUHF)
            CALL CMPENG(ICORE(I0),MAXCOR,60,0,ECORR,ENERGY(2,1),
     &                  ENERGY(2,2),IUHF,0)
            WRITE(6,4001)ENERGY(2,1)
4001        FORMAT(T3,'E3D = ',F14.10,' a.u.')
4002        FORMAT(T3,'E4D = ',F14.10,' a.u.')
4003        FORMAT(T3,'E4Q = ',F14.10,' a.u.')
C
C  FOR GRADIENT CALCULATION MOVE T2 TO THE LAMBDA LISTS
C
            IF (GRAD) THEN
               CALL INIT2(IUHF)
               CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,.FALSE.)
            END IF
            WRITE(6,4002)ENERGY(3,1)
            CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
            CALL GENINT(ICORE(I0),MAXCOR,IUHF,1)
            CALL FEACONT(ICORE(I0),MAXCOR,IUHF)
            CALL FMICONT(ICORE(I0),MAXCOR,IUHF)
            CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,1)
            CALL NEWT2(ICORE(I0),MAXCOR,IUHF)
            CALL CMPENG(ICORE(I0),MAXCOR,60,0,ECORR,ENERGY(4,1),
     &                  ENERGY(4,2),IUHF,0)
            WRITE(6,4003)ENERGY(4,1)
            CALL E4S(ICORE(I0),MAXCOR,IUHF,ENERGY(5,1))
            IF (TRIPNI) CALL TRPS(ICORE(I0),MAXCOR,IUHF,ENERGY(6,1))
            IF (IUHF.EQ.0) CALL RESET(ICORE(I0),MAXCOR,IUHF)

c        END IF (NONHF)
         END IF

         IF (IFLAGS(35).NE.0) CALL ACES_AUXCACHE_FLUSH
         CALL MBPTOUT

c     ELSE IF (CC) THEN
      ELSE

      do 2341 iref=1,nref
cjp parallelization
      if (mod(iref-1,nproc).ne.myproc-1) goto 2341
c
      if (nref.gt.nproc) call reopenmo
      CALL INMBPT(ICORE(I0),MAXCOR,IUHF)
C
C CC CALCULATION LOGIC BEGINS HERE.
C
C SG 7/24/98
C Try to read in a T guess.  If we can, set SING1 to be true.
C
      CALL DRDTGSS(ICORE(I0),MAXCOR/IINTFP,IUHF,'TGUESS  ',0,READT)
C
C Note that we need to test to see whether we are doing a CCD
C calculation before we turn the SING1 on. This was not tested
C originaly and CCD optimizations were failing. 01/2001, Ajith Perera
C
      IF ((IFLAGS(2).NE.5.AND.IFLAGS(2).NE.8).AND.READT) SING1 = .TRUE.
C
      CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
C
C GENERATE W LISTS FOR INTERMEDIATES IF THIS IS NOT MBPT(3) OR LCCD.
C
C      IF(METHOD.NE.2.AND.METHOD.NE.5)THEN
      CALL SETLST(ICORE(I0),MAXCOR,IUHF)
C      ENDIF
C**********************************************************************
C
C START OF CC LOOP
C
      RLECYC = 0
      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(ICORE(I0),MAXCOR,IUHF,RLECYC,.FALSE.)
      END IF
      IF (IFLAGS(21).EQ.1) THEN
cjp diis probably should be iref specific
         CALL DIISLST(1,IUHF,METHOD.GE.6.AND.METHOD.NE.8)
      END IF
cjp@@@ here we have a problem: if we allow some stuff to be in core,
cjp       if the flushcache will really synchronize it into  file!
      CALL INCOR(I0,IUHF)
      MAXCOR = ICRSIZ
      if (nref.gt.nproc) then
         call aces_cache_flush
         call storemoio
      end if

2341  continue
cjp end of iref loop 2

      idiisstop = 0
cjp
cjp allow for ncycle change during iterations, avoid do loop
cjp - change of loop limits inside loop might be non-portable
c
c      DO 55 ICYCLE=1,NCYCLE
      icycle = 1
555   continue
      RLECYC = RLECYC + 1
C      IF (ICYCLE.EQ.1) CALL UPDVC(ICORE(I0),ECORR,1,NCYCLE,IUHF)
C
C FILL T2 INCREMENTS WITH <IJ||AB> INTEGRALS
C
C AND ZERO THE T1 INCREMENTS (ONLY CCSD AND QCISD)
C
CJDW 10/24/95. Last argument was removed from INITIN, so remove it here.
C      CALL INITIN(ICORE(I0),MAXCOR,IUHF,IFLAGS(93).EQ.2)

      do 3412 iref = 1, nref
cjp parallelization
      if (mod(iref-1,nproc).ne.myproc-1) goto 3412
c
      if (nref.gt.nproc) call reopenmo
      CALL INITIN(ICORE(I0),MAXCOR,IUHF)
      IF (METHOD.GE.6.AND.METHOD.NE.8) THEN
         CALL INITSN(ICORE(I0),MAXCOR,IUHF)
      END IF
C
C GENERATE W AND F INTERMEDIATES.
C
      IF (LINCC) THEN
         INTTYP = 3
      ELSE
         INTTYP = 2
      END IF
      CALL GENINT(ICORE(I0),MAXCOR,IUHF,INTTYP)
C
C COMPUTE F INTERMEDIATE CONTRIBUTION TO T2 INCREMENT.
C
      CALL FEACONT(ICORE(I0),MAXCOR,IUHF)
      CALL FMICONT(ICORE(I0),MAXCOR,IUHF)
      IF ((METHOD.GT.9.AND.SING1).OR.(METHOD.EQ.6.AND.SING1)) THEN
         CALL FMECONT(ICORE(I0),MAXCOR,IUHF,1)
         IF (IUHF.NE.0) CALL FMECONT(ICORE(I0),MAXCOR,IUHF,2)
         CALL T1INT2A(ICORE(I0),MAXCOR,IUHF)
         CALL T1INT2B(ICORE(I0),MAXCOR,IUHF)
         CALL T1INT1(ICORE(I0),MAXCOR,IUHF,1)
         IF (IUHF.NE.0) CALL T1INT1(ICORE(I0),MAXCOR,IUHF,2)
      END IF
      IF (UCC) THEN
         CALL T1INT2A(ICORE(I0),MAXCOR,IUHF)
         CALL T1INT2B(ICORE(I0),MAXCOR,IUHF)
      END IF
C
C DO W INTERMEDIATE CONTRIBUTION TO T2 EQUATION
      CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)
C
C     DO TdegerWT CONTRIBUTION FOR UCC AND XCC
      IF (UCC) CALL TTWT(ICORE(I0),MAXCOR,IUHF)
C
C NOW DO T3->T2 AND T3->T1 FOR ITERATIVE TRIPLE METHODS
      IF (TRIPIT) THEN
         CALL TRPINT(ICORE(I0),MAXCOR,IUHF)
         CALL TRPS(ICORE(I0),MAXCOR,IUHF,ZJUNK)
      END IF
      IF (OSFLAG) CALL OS(ICORE(I0),MAXCOR,IUHF,W)
C
C COMPUTE THE CONTRIBUTION OF DOUBLES TO T1 AND
C DENOMINATOR WEIGHT THE NEW INCREMENTS
C (ONLY CCSD AND QCISD METHODS)
C
      IF (METHOD.GE.9.OR.METHOD.EQ.6.OR.METHOD.EQ.7) THEN
cjp set flag to perform only first part of processing in e4s, e4seng
         if (isbwcc) then
            ibwpass=1
         else
            ibwpass=0
         end if
         CALL E4S(ICORE(I0),MAXCOR,IUHF,EDUMMY)
      END IF
C
C DENOMINATOR WEIGHT T2 INCREMENTS TO FORM NEW T2.
C
cjp NEWT2 uses common bwcc and conditionals
cjp set flag to perform only first part of processing
      if (isbwcc) then
         ibwpass=1
      else
         ibwpass=0
      end if
      CALL NEWT2(ICORE(I0),MAXCOR,IUHF)
cjp
      if (nref.gt.nproc) then
         call aces_cache_flush
         call storemoio
      end if

3412  continue

cjp
      if (isbwcc) then
cjp communicate columns of heff between processors
cjp in parallel calculation
cjp diagonalize heff
         call heffdiag(iuhf,icycle)
         if (myproc.eq.1) then
            write(6,*)'@VCC: MR-BW-CC CORRELATION ENERGY = ',epsilon0
            write(6,*)
         end if
      end if
cjp

cjp store correctly in energy array, ecorr variable
      ixcycle = -1
      do 4123 iref = 1, nref
cjp parallelization
      if (mod(iref-1,nproc).ne.myproc-1) goto 4123
      if (nref.gt.nproc) call reopenmo

      if (isbwcc) then
cjp setup BW shift of denominator
         if (ihubaccorr.lt.2) then
cjp ordinary BWCC value
            ecorrbw=epsilon0-heff(iref,iref)
         else
cjp FOR a posteriori size extenzivity correction make it like in RSPT-CC
            ecorrbw=0.d0
         end if
         ecorrbw0=epsilon0-heff(iref,iref)
cjp homotopical transition between BW and RS
         hfakt=1.0
         if (ihomotop.gt.0) then
            if (icycle.gt.ihomotop) then
cjp first such iteration ... switch off diis
               if (icycle.eq.ihomotop+1) idiisstop = 1
               hfakt=lambdahomotop**(icycle-ihomotop)
cjp restart diis eventually
               if (hfakt.lt.diishonset.and.
     &             idiisstop.eq.1.and.
     &             icycle.gt.ihomotop+1) then
cjp postpone this in order to call diislst for other references, too
cjp but set a flag that in this cycle, dodiis0 should already be called
                  ixcycle=icycle
                  if (iref+nproc.gt.nref) idiisstop=0
cjp remember, we are still inside iref loop
                  IF (IFLAGS(21).EQ.1) THEN
                     write(6,*)'DIIS RESTART',icycle,hfakt,diishonset
                     CALL DIISLST(1,IUHF,METHOD.GE.6.AND.METHOD.NE.8)
                  END IF
               end if
            end if
            ecorrbw0=epsilon0-heff(iref,iref)
            ecorrbw=hfakt*ecorrbw0
            if (icycle.gt.ihomotop) then
               write(6,*) 'homotopic solution: ',
     &                'cycle,hfakt,scaled ecorrbw=',icycle,hfakt,ecorrbw
            end if
         end if
cjp here we must finish the work of e4seng, newt2
         ibwpass=2
         call newt2(icore(i0),maxcor,iuhf)
cjp NOTE: PUTLST(ICORE(ID),1,1,1,3,90) contains elements needed for  eq 4.28
cjp which will be destroyed and replaced by final t1 amplitudes in e4s call
cjp thus newt2 is here called first, in contrast to original code
cjp@@@!!!note they correspond to NEW t1 amplitudes!
         call e4s(icore(i0),maxcor,iuhf,edummy)
c     end if (isbwcc)
      end if
C
C WRITE OUT RMS AND MAX DIFFERENCES FOR NEW T2 AND PUT NEW ERROR VECTOR
      if (isbwcc) then
         call drtsts(icore(i0),maxcor,icycle,iuhf,
     &               ibwconvg(iref),icontl,sing1,0,'t')
      else
         CALL DRTSTS(ICORE(I0),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &               SING1,0,'T')
      end if
cjp this computes already new diagonal heff elements to be used in next iter
cjp NOTE this overwrites heff
      if (isbwcc) then
cjp in energy array store the true mr-cc energy!
         call cmpeng(icore(i0),maxcor,60,2,ecorr,
     &               edummy1,edummy2,iuhf,1)
         energy(icycle+1,1)=epsilon0
         energy(icycle+1,2)=epsilon0
      else
         CALL CMPENG(ICORE(I0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,1)
      end if
      IF (OSFLAG) THEN
         WRITE(6,81)ENERGY(ICYCLE+1,1)+W
 81      FORMAT(T3,' The total correlation energy for the singlet ',
     &          'is ',F15.12,' a.u.')
         WRITE(6,82)ENERGY(ICYCLE+1,1)-W
 82      FORMAT(T3,' The total correlation energy for the triplet ',
     &          'is ',F15.12,' a.u.')
         ENERGY(ICYCLE+1,1)=ENERGY(ICYCLE+1,1)+W
         ENERGY(ICYCLE+1,2)=ENERGY(ICYCLE+1,2)+W
      END IF
      ELAST = ENERGY(ICYCLE+1,2)
C
C IF ICONVG IS EQUAL ZERO, CONVERGENCE HAS BEEN ACHIEVED, EXIT
C VIA FINISH
C
cjp switch off diis for a posteriori correction
      IF (IFLAGS(21).EQ.1.and.ihubaccorr.lt.2) THEN
         if (.not.isbwcc) then
            CALL DODIIS0(ICORE(I0),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                   ICONVG,ICONTL,SING1,44,61,90,0,90,2,70,
     &                   '     ')
         else
cjp bwcc case
cjp NOTE: strictly taken, we do illegal thing by calling diis INDEPENDENTLY
cjp for each iref,
cjp we should better build a common vector for all irefs and do ONE diis
cjp however, in the cases test up to now it did not matter too much, perhaps
cjp in last iterations with strict conv. criterion is wastes 1-2 iterations
cjp
            if (idiisstop.eq.0 .or. icycle.eq.ixcycle) then
               call dodiis0(icore(i0),maxcor/iintfp,iuhf,1,icycle,
     &                      ibwconvg(iref),icontl,sing1,44,61,90,0,90,2,
     &                      70,  '     ')
            end if
         end if
      END IF
c
cjp NOTE this iref loop has been split in order to run the amplitude
cjp update consistently for all references
cjp
      if (nref.gt.nproc) then
         call aces_cache_flush
         call storemoio
      end if

4123  continue
c
cjp TEST OF CONVERGENCE SHOULD BE OUTSIDE THE IREF LOOP
      if (isbwcc) then
#ifdef MPI
cjp communicate preliminary results for all references
         if (nproc.gt.1) then
            do iref = 1, nref
               call MPI_BCAST(ibwconvg(iref),1,MPI_INTEGER,
     &                        mod(iref-1,nproc),MPI_COMM_WORLD,ierror)
            end do
         end if
#endif
         if (ihomotop.gt.0.and.hfakt.gt.hfaktmax) then
cjp force further iterations if homotopy is not completed yet
            iconvg=1
         else
cjp test now on all processors
            iconvg=0
         end if
cjp iconvg zero, if all of them zero for all configs
         do iref = 1, nref
            if (ibwconvg(iref).ne.0) iconvg=1
         end do
      end if
cjp
cjp
cjp here is a place to exit from mr-bwcc loop
      if (isbwcc.and.iconvg.eq.0) then
         write(6,*)'@VCC: the MR-BW-CC calculation has converged'
      end if
cjp
cjp implementation of the Hubac size extenzivity correction for BWCC
      if (isbwcc.and.ihubaccorr.eq.1.and.iconvg.eq.0.or.
     &    isbwcc.and.ihubaccorr.gt.1) then
cjp pass one ... switch off diis for the future, compute another iteration in
cjp which RSPT denominators will be used ... pass 2
cjp in pass 3, the new amplitudes obtained using RSPT denominators will be used
cjp to obtain new Heff and the final energy by its diagonalization
         if (ihubaccorr.eq.3) then
cjp finish!
            iconvg=0
         end if
         if (ihubaccorr.eq.2) then
            iconvg=1
            ihubaccorr=3
         end if
         if (ihubaccorr.eq.1) then
            write(6,*) 'Computing the size-extenzivity correction'
            iconvg=1
            ihubaccorr=2
         end if
      end if
cjp
cjp start homotopy if convergence has been reached
cjp
      if (isbwcc.and.ihomotop.gt.0.and.iconvg.eq.0.and.
     &    icycle.lt.ihomotop) then
         ihomotop=icycle
         iconvg=1
      end if
cjp
      do 4321 iref=1,nref
cjp parallelization
      if (mod(iref-1,nproc).ne.myproc-1) goto 4321
      if (nref.gt.nproc) call reopenmo
C
CJDW 10/24/95.
C    Do not update if we have converged. This allows us to maintain a
C    correspondence between the amplitudes and existing Hbar elements.
C    This is important for MN's finite-order EOM stuff.
C
cjp for mrbwcc, we also do not overwrite the amplitudes, from which the current
cjp energy was computed
cjp use the all-iref iconvg here!
      IF (ICONVG.NE.0) THEN
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      END IF
CJDW END
cjp when convergence was achieved
      IF (ICONVG.EQ.0) THEN
         if (isbwcc) then
cjp the cmpeng call is probably obsolete here
cjp  call cmpeng(icore(i0),maxcor,43,0,ecorr,edummy1,edummy2,iuhf,1)
            energy(icycle+1,1)=epsilon0
            energy(icycle+1,2)=epsilon0
         else
            CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,ENERGY(ICYCLE+1,1),
     &                  ENERGY(ICYCLE+1,2),IUHF,1)
         end if
c& new lines
         IF (OSFLAG) THEN
            WRITE(6,81)ENERGY(ICYCLE+1,1)+W
            WRITE(6,82)ENERGY(ICYCLE+1,1)-W
            ENERGY(ICYCLE+1,1)=ENERGY(ICYCLE+1,1)+W
            ENERGY(ICYCLE+1,2)=ENERGY(ICYCLE+1,2)+W
         END IF
c& new lines end
      END IF
cjp        ENDIF
C SG 7/21/98 Check if we are ready for a Brueckner step
      IF (BRUECK .AND. (ICYCLE.GT.1)) THEN
         CALL CHKBKSTP(ICORE(I0), MAXCOR/IINTFP, IUHF, ICONVG)
      END IF
C
      IF (ICONVG.EQ.0) THEN
C
C  CALCULATE THE UCC AND XCC ENERGY CONTRIBUTION
C   !!!! NOTE THAT THIS ROUTINE REWRITES WMBEJ LISTS 54-59 !!!!
C   !!!! NOTE THAT THIS ROUTINE REWRITES T2INC LISTS 41-46 !!!!
         IF (UCC) CALL EUCC(ICORE(I0),MAXCOR,IUHF,ENERGY(ICYCLE+1,1))
C
CJDW 3/20/96. Print final amplitudes.
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
C
C  PUT CONVERGED AMPLITUDES AS INITIAL GUESS FOR SOLVING THE
C  LAMBDA EQUATIONS ON THE GAMLAM FILE
         WRITE(6,4000)
 4000    FORMAT(T3,' The CC iterations have converged.',/)
         IF ((GRAD.AND..NOT.UCC).OR.
     &       (IFLAGS(87).EQ.3 .AND. (TRIPIT.OR.TRIPNI1) ) ) THEN
            CALL INIT2(IUHF)
            CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,100,SING1)
         END IF
         IF (IUHF.NE.0) CALL S2PROJ(ICORE(I0),MAXCOR,IUHF,SING1)
C
C-----------------------------------------------------------------------
C     In EOM triples calculations, RESET is sometimes called by RESRNG.
C     Determine if we need to call RESET here or not.
C-----------------------------------------------------------------------
         DORESET = IFLAGS(87).NE.3 .OR.
     &            (IFLAGS(87).EQ.3 .AND. IFLAGS(2).NE.16 .AND.
     &                                   IFLAGS(2).NE.18 .AND.
     &                                   IFLAGS(2).NE.33 .AND.
     &                                   IFLAGS(2).NE.34 .AND.
     &           .NOT. (IFLAGS(2).EQ.22 .AND. IFLAGS2(124).GE.5))
         IF (DORESET) THEN
cjp this should be called for all iref
cjp leave it here in the iref loop 5
            IF(IUHF.EQ.0.AND..NOT.UCC)CALL RESET(ICORE(I0),MAXCOR,IUHF)
         END IF
C
         IF (IFLAGS(87).EQ.3.AND.(TRIPIT.OR.TRIPNI)) THEN
            CALL RESRNG(ICORE(I0),MAXCOR,IUHF)
         END IF
C-----------------------------------------------------------------------
C
cjp this has sense only to be called once!
         if (.not.isbwcc .or. iref+nproc.gt.nref) then
            CALL FINISH(ICYCLE+1)
         end if
C SG 7/21/98
         CALL DDMPTGSS(ICORE(I0), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')
C
         IF (BRUECK) CALL BRUECKIT(ICORE(I0),MAXCOR,IUHF,BRKCNV)
         IF (TRIPNI .AND. (.NOT.BRUECK .OR. BRKCNV)) THEN
            CALL TRPS(ICORE(I0),MAXCOR,IUHF,EDUMMY)
         END IF
         IF (IFLAGS(35).NE.0) CALL ACES_AUXCACHE_FLUSH
C
C-----------------------------------------------------------------------
C     Close the T3 file. Delete it except in certain noniterative meth-
C     ods and in EOMEE triples methods.
C-----------------------------------------------------------------------
C
         IF (T3STOR) THEN
            IF ((METACT.EQ.29).OR.(METACT.EQ.30)) METHOD = METACT
            IF ((METHOD.EQ.12).OR.(METHOD.GE.26.AND.METHOD.LE.31).OR.
     &          (IFLAGS(87).EQ.3.AND.(TRIPIT.OR.TRIPNI))) THEN
               CALL AUXIOO(.FALSE.)
            ELSE
               CALL AUXIOO(.TRUE.)
            END IF
         END IF
C
C-----------------------------------------------------------------------
CJDW 10/24/95. Write energy and coordinates in POLYRATE jobs.
C-----------------------------------------------------------------------
         IF (IFLAGS2(113).GT.0) THEN
            IONE=1
            CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOMS)
            CALL POLYPRT0(NATOMS,IINTFP,IFLAGS(1),
     &                    .TRUE.,.TRUE.,.FALSE.,.FALSE.)
         END IF
C-----------------------------------------------------------------------
         call aces_fin
cjp go through whole iref loop
         if (.not.isbwcc .or. iref+nproc.gt.nref) then
cjp ... for eneronly geometry optimizations - either in each program
cjp     gfname must know that MOINTS are MOINTS.n
cjp     or simpler and at the moment sufficient, MOINTS will be deleted
            if (isbwcc) then
               istat = ishell('rm -f MOINTS.* MOABCD.*')
               write(6,*) '@VCC-I: files MOINTS.* MOABCD.* removed'
            end if
            WRITE(6,1020)
            call bwparend
            call aces_exit(0)
         end if

c     END IF (ICONVG.EQ.0)
      END IF

C DO THE RLE EXTRAPOLATION AND GET THE UPDATED RLE ENERGY.
C
cjp no EXTRAPOLATION for the a posteriori correction
      IF (MOD(IFLAGS(21),2).EQ.0 .and. ihubaccorr.lt.2) THEN
         CALL DRRLE(ICORE(I0),MAXCOR,IUHF,RLECYC,.FALSE.)
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      END IF
C
C PRINT AMPLITUDES EVERY IEVERY CYCLES.
C
cjp more default print for debug
cjp print the true converged amplitudes if correction is to be calculated
cjp since they will be overwritten by amplitudes obtained using the
cjp modified denominator
      IF (iconvg.ne.0 .and. ihubaccorr.ge.2 .or.
     &    isbwcc.and.(bwgossip.or.ihomotop.gt.0).and.iconvg.ne.0.or.
     &    IFLAGS(1).GE.10 .and. iconvg.ne.0) THEN
         IF (ihubaccorr.ge.2 .or. MOD(ICYCLE,IEVERY).EQ.0) THEN
            if (ihubaccorr.eq.2) then
               write(6,*) ' Ordinary BWCC converged amplitudes'
            end if
            if (ihubaccorr.eq.3) write(6,*) ' Corrected BWCC amplitudes'
            CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
         END IF
      END IF
cjp
      if (nref.gt.nproc) then
         call aces_cache_flush
         call storemoio
      end if

4321  continue
cjp
cjp55     CONTINUE
      icycle=icycle+1
      if (icycle.le.ncycle) goto 555
cjp
      WRITE(*,*)
     &         '@VCC: The Coupled-Cluster equations did not converge!!!'
cYAU - If we are here, did not the CC equations fail to converge?
      call aces_exit(1)
C***********************************************************************
      do iref = 1, nref
         if (mod(iref-1,nproc).eq.myproc-1) then
            if (nref.gt.nproc) call reopenmo
            CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
            if (nref.gt.nproc) then
               call aces_cache_flush
               call storemoio
            end if
         end if
      end do

c     END IF (.NOT.CC)
      END IF

CJDW 10/24/95. Write energy and coordinates in POLYRATE jobs.
      IF (IFLAGS2(113).GT.0) THEN
         IONE=1
         CALL GETREC(20,'JOBARC','NATOMS  ',IONE,NATOMS)
         CALL POLYPRT0(NATOMS,IINTFP,IFLAGS(1),
     &                 .TRUE.,.TRUE.,.FALSE.,.FALSE.)
      END IF

C IN THE CASE OF OUT OF CORE ALGORITHM SAVE AMPLITUDE AND INTERMEDIATE
C LISTS IF NECESSARY
      call aces_fin

 1020 FORMAT(/,77('-'),/,32X,'Exiting xvcc',/,77('-'),/)
      WRITE(6,1020)
      call bwparend

      END

