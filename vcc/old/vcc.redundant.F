
C VECTORIZED CC/MBPT PROGRAM
C
C WRITTEN BY J.F. STANTON, J. GAUSS, J.D. WATTS AND W.J. LAUDERDALE
C TD-CC EXTENSION ADDED BY P.G. SZALAY (ROUTINES WITH 'OS')

      PROGRAM VCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISTSZ,RLECYC
      INTEGER POP(8,2),VRT(8,2)
      logical bPrint
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
#include "bwcc.com"
      integer ncycle
      common /ncycle/ncycle
C
      EQUIVALENCE (IFLAGS(2),METHOD)
C
 1000 FORMAT(/,77('-'),/,32X,'Entering xvcc',/,77('-'),/)
      WRITE(6,1000)
C
      call bwread(.true.,iuhf)
      if (isbwcc) call bwcc
C
      iref = 1
      ibwpass = 0

      CALL CRAPSI(ICORE,IUHF,0)
      bPrint = iflags(1).ge.10

      CALL ZERO(ENERGY,1000)
      CALL IZERO(IXTRLE,500)
      CALL SETMET
      CALL SETLOG(RESTART)
      ICONVG=1
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

C INITIALIZE DOUBLES LISTS
      IF (METHOD.GT.1) CALL INMBPT(ICORE(I0),MAXCOR,IUHF)

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
C
      IF (CC) THEN
         NCYCLE = IFLAGS(7)
         NKEEP  = IFLAGS(12)
         ICONTL = IFLAGS(4)
      END IF
      IF (bPrint) THEN
         WRITE(6,801)
  801    FORMAT(T3,' Initial T amplitudes: ')
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
         WRITE(6,800)
  800    FORMAT(T3,' Correlation energies computed from initial ',
     &          'T amplitudes: ')
      END IF
      CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,
     &            ENERGY(1,1),ENERGY(1,2),IUHF,1)
      ELAST = ENERGY(1,2)

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
            stop
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
C FOURTH-ORDER LOGIC.  COMPUTE D3 AND D4 ON FIRST PASS, THEN GET QUADS.
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
C Start the T2 lists by adding <ab||ij>, This is actually a 
C contribution to T2.
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
         CALL DIISLST(1,IUHF,METHOD.GE.6.AND.METHOD.NE.8)
      END IF
      CALL INCOR(I0,IUHF)
      MAXCOR = ICRSIZ

      icycle = 1
555   continue
C
      RLECYC = RLECYC + 1
C
C IF (ICYCLE.EQ.1) CALL UPDVC(ICORE(I0),ECORR,1,NCYCLE,IUHF)
C FILL T2 INCREMENTS WITH <IJ||AB> INTEGRALS
C AND ZERO THE T1 INCREMENTS (ONLY CCSD AND QCISD)
C JDW 10/24/95. Last argument was removed from INITIN, so remove it here.
C CALL INITIN(ICORE(I0),MAXCOR,IUHF,IFLAGS(93).EQ.2)
C
C Start the T2 lists by adding <ab||ij>, This is actually a
C contribution to T2.
C
      CALL INITIN(ICORE(I0),MAXCOR,IUHF)
C
      IF (METHOD.GE.6.AND.METHOD.NE.8) THEN
C
C Zero out the T1 lists for all cycles.
C
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
C 
C Generate W(mbej), W(mnij) W(mbej, F(ea), F(ij) and F(me)
C intermediates.  

      CALL GENINT(ICORE(I0),MAXCOR,IUHF,INTTYP)
C
C COMPUTE F INTERMEDIATE CONTRIBUTION TO T2 INCREMENT.
C
C Evaluate P(ab)Sum_e T(ij,ae){F(b,e)-1/2 sum_m T(m,b)F(m,e)}
C contribution to T2 and Sum_e T(i,e)F(a,e) contribution to T1.
C
      CALL FEACONT(ICORE(I0),MAXCOR,IUHF)
C
C Evaluate -P(ij)Sum_e T(im,ab){F(mj)-1/2 sum_m T(j,e)F(m,e)}
C contribution to T2 and Sum_e T(m,a)F(m,i) contribution to T1.
C
      CALL FMICONT(ICORE(I0),MAXCOR,IUHF)
C
      IF ((METHOD.GT.9.AND.SING1).OR.(METHOD.EQ.6.AND.SING1)) THEN
C
C Evaluate Sum_e T(im,ae)F(m,e) contribution to T1. 
C
         CALL FMECONT(ICORE(I0),MAXCOR,IUHF,1)
C
         IF (IUHF.NE.0) CALL FMECONT(ICORE(I0),MAXCOR,IUHF,2)
C 
C Evaluate -P(ab)Sum_m T(m,a)<mb||ij> contribution to T2.
C 
         CALL T1INT2A(ICORE(I0),MAXCOR,IUHF)
C
C Evaluate +P(ij)Sum_m T(i,e)<ab||ej> contribution to T2.
C
         CALL T1INT2B(ICORE(I0),MAXCOR,IUHF)
C
C Evaluate -Sum_nf T(n,f)<na||if> contribution to T1.
C
         CALL T1INT1(ICORE(I0),MAXCOR,IUHF,1)
         IF (IUHF.NE.0) CALL T1INT1(ICORE(I0),MAXCOR,IUHF,2)
C
      END IF
C 
      IF (UCC) THEN
         CALL T1INT2A(ICORE(I0),MAXCOR,IUHF)
         CALL T1INT2B(ICORE(I0),MAXCOR,IUHF)
      END IF
C
C DO W INTERMEDIATE CONTRIBUTION TO T2 EQUATION
C 
C Evaluate the following three contributions to T2:
C
C   -P(ij)P(ab) Sum_me T(i,e)T(m,a) <mb||ej> (T12INT2)
C   +P(ij)P(ab) Sum_me T(i,e)T(m,a) W(mb,ej) (DRRNG)
C   +1/2 Sum_mn Tau(mn,ab)W(mn,ij) + 
C    1/2 Sum_ef Tau(ij,ef)W(ab,ef)           (DRLAD)
C
      CALL DRE3EN(ICORE(I0),MAXCOR,IUHF,0)
C
C DO TdegerWT CONTRIBUTION FOR UCC AND XCC
C 
      IF (UCC) CALL TTWT(ICORE(I0),MAXCOR,IUHF)
C
C NOW DO T3->T2 AND T3->T1 FOR ITERATIVE TRIPLE METHODS
C
      IF (TRIPIT) THEN
         CALL TRPINT(ICORE(I0),MAXCOR,IUHF)
         CALL TRPS(ICORE(I0),MAXCOR,IUHF,ZJUNK)
      END IF
C 
      IF (OSFLAG) CALL OS(ICORE(I0),MAXCOR,IUHF,W)
C
C COMPUTE THE CONTRIBUTION OF DOUBLES TO T1 AND
C DENOMINATOR WEIGHT THE NEW INCREMENTS
C (ONLY CCSD AND QCISD METHODS)
C 
      IF (METHOD.GE.9.OR.METHOD.EQ.6.OR.METHOD.EQ.7) THEN
C
C Evaluate the following three contributions to T1:
C
C   -1/2 Sum_mef T(im,ef)<ma||ef>  (T2T1AA1, T2T1AB1)
C   -1/2 Sum_men T(nm,ei)<nm||ei>  (T2T1AA2, T2T1AB2)
C
C Also, Do the T1 = T1/{f(i,i) - f(a,a)} to get a new T1
C
         CALL E4S(ICORE(I0),MAXCOR,IUHF,EDUMMY)
      END IF
C
C DENOMINATOR WEIGHT T2 INCREMENTS TO FORM NEW T2.
C 
C Do the T2 = T2/{f(i,i) + f(j,j) - f(a,a) - f(b,b)}
C to get a new T2.
C
      CALL NEWT2(ICORE(I0),MAXCOR,IUHF)
C
C WRITE OUT RMS AND MAX DIFFERENCES FOR NEW T2 AND PUT 
C NEW ERROR VECTOR
C
      CALL DRTSTS(ICORE(I0),MAXCOR,ICYCLE,IUHF,ICONVG,ICONTL,
     &            SING1,0,'T')
      CALL CMPENG(ICORE(I0),MAXCOR,60,2,ECORR,ENERGY(ICYCLE+1,1),
     &            ENERGY(ICYCLE+1,2),IUHF,1)
C
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
      IF (IFLAGS(21).EQ.1) THEN
         CALL DODIIS0(ICORE(I0),MAXCOR/IINTFP,IUHF,1,ICYCLE,
     &                ICONVG,ICONTL,SING1,44,61,90,0,90,2,70,
     &                '     ')
      END IF

CJDW 10/24/95.
C    Do not update if we have converged. This allows us to maintain a
C    correspondence between the amplitudes and existing Hbar elements.
C    This is important for MN's finite-order EOM stuff.
C
      IF (ICONVG.NE.0) THEN
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      END IF
CJDW END

      IF (ICONVG.EQ.0) THEN
         CALL CMPENG(ICORE(I0),MAXCOR,43,0,ECORR,ENERGY(ICYCLE+1,1),
     &               ENERGY(ICYCLE+1,2),IUHF,1)
         IF (OSFLAG) THEN
            WRITE(6,81)ENERGY(ICYCLE+1,1)+W
            WRITE(6,82)ENERGY(ICYCLE+1,1)-W
            ENERGY(ICYCLE+1,1)=ENERGY(ICYCLE+1,1)+W
            ENERGY(ICYCLE+1,2)=ENERGY(ICYCLE+1,2)+W
         END IF
      END IF

C SG 7/21/98 Check if we are ready for a Brueckner step
      IF (BRUECK .AND. (ICYCLE.GT.1)) THEN
         CALL CHKBKSTP(ICORE(I0), MAXCOR/IINTFP, IUHF, ICONVG)
      END IF

      IF (ICONVG.EQ.0) THEN

C  CALCULATE THE UCC AND XCC ENERGY CONTRIBUTION
C   !!!! NOTE THAT THIS ROUTINE REWRITES WMBEJ LISTS 54-59 !!!!
C   !!!! NOTE THAT THIS ROUTINE REWRITES T2INC LISTS 41-46 !!!!
         IF (UCC) CALL EUCC(ICORE(I0),MAXCOR,IUHF,ENERGY(ICYCLE+1,1))

CJDW 3/20/96. Print final amplitudes.
         CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')

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
            IF(IUHF.EQ.0.AND..NOT.UCC)CALL RESET(ICORE(I0),MAXCOR,IUHF)
         END IF
C
         IF (IFLAGS(87).EQ.3.AND.(TRIPIT.OR.TRIPNI)) THEN
            CALL RESRNG(ICORE(I0),MAXCOR,IUHF)
         END IF
C-----------------------------------------------------------------------

         CALL FINISH(ICYCLE+1)

C SG 7/21/98
         CALL DDMPTGSS(ICORE(I0), MAXCOR/IINTFP, IUHF, 0, 'TGUESS  ')

         IF (BRUECK) CALL BRUECKIT(ICORE(I0),MAXCOR,IUHF,BRKCNV)
         IF (TRIPNI .AND. (.NOT.BRUECK .OR. BRKCNV)) THEN
            CALL TRPS(ICORE(I0),MAXCOR,IUHF,EDUMMY)
         END IF
         IF (IFLAGS(35).NE.0) CALL ACES_AUXCACHE_FLUSH

C-----------------------------------------------------------------------
C     Close the T3 file. Delete it except in certain noniterative meth-
C     ods and in EOMEE triples methods.
C-----------------------------------------------------------------------
         IF (T3STOR) THEN
            IF ((METACT.EQ.29).OR.(METACT.EQ.30)) METHOD = METACT
            IF ((METHOD.EQ.12).OR.(METHOD.GE.26.AND.METHOD.LE.31).OR.
     &          (IFLAGS(87).EQ.3.AND.(TRIPIT.OR.TRIPNI))) THEN
               CALL AUXIOO(.FALSE.)
            ELSE
               CALL AUXIOO(.TRUE.)
            END IF
         END IF

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
         WRITE(6,1020)
         stop

c     END IF (ICONVG.EQ.0)
      END IF

C DO THE RLE EXTRAPOLATION AND GET THE UPDATED RLE ENERGY.
C
      IF (MOD(IFLAGS(21),2).EQ.0) THEN
         CALL DRRLE(ICORE(I0),MAXCOR,IUHF,RLECYC,.FALSE.)
         CALL DRMOVE(ICORE(I0),MAXCOR,IUHF,0,SING1)
         CALL RNABIJ(ICORE(I0),MAXCOR,IUHF,'T')
      END IF
C
C PRINT AMPLITUDES EVERY IEVERY CYCLES.
C
      IF (bPrint .and. iconvg.ne.0) THEN
         IF (MOD(ICYCLE,IEVERY).EQ.0) THEN
            CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
         END IF
      END IF

      icycle=icycle+1
      if (icycle.le.ncycle) goto 555

      WRITE(*,*) '@VCC: The Coupled-Cluster equations did not converge!'
      CALL AMPSUM(ICORE(I0),MAXCOR,IUHF,0,SING1,'T')
      call aces_exit(1)

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

      END

