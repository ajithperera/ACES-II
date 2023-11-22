#include "flags.h"
      SUBROUTINE DOPARTEOM(IUHF, SCR, MAXCOR)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION ICONV(100), IREORD(100)
      DIMENSION SCR(MAXCOR)
      LOGICAL ESPROP
      LOGICAL LAMTHERE, TRIPLET, DOBLOCK, LOCKLEFT, SHIFTROOT
      LOGICAL MBPT2, CC, SS, SD, DS, DD, CCD, LCCD,LCCSD
      LOGICAL RCCD,DRCCD
      LOGICAL CORE_PROJECT,CORE_SEARCH,DOUBLE_EES
      INTEGER DIRPRD
      INTEGER END,BGN,BGN_IRP,END_IRP
      DOUBLE PRECISION EIGVAL_T
      CHARACTER*5 HAND(2)
      CHARACTER*1 NATURE(100,8)
C
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/CALCINFO/NROOT(8)
      COMMON/ROOTS/EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),
     &             BGN(100,8),BGN_IRP(100,8),END(100,8),
     &             END_IRP(100,8),NATURE
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC 
      COMMON/PROPGRAD/ESPROP,IOPTROOT,IOPTSYM
      COMMON/EIGPROB/ISIDE
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD
      COMMON/LISTPROJ/LISTH0, ICOLPR1, ICOLPR2
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON/DRVHBAR/SS, SD, DS, DD
      COMMON/TDALIST/LISTETDA, LISTVTDA
      COMMON/LAMSTATE/LAMTHERE
      COMMON/SPINSTATE/TRIPLET
      COMMON/SHFTROOT/SHIFTROOT
      COMMON/STATSYM/IRREPX
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA ONE /1.0D0/
      DATA FACTEV /27.2113957D0/
      DATA HAND /'right',' left'/
C
      LISTRIGHT = 495
      LISTLEFT  = 496
      LISTE     = 497
      LISTPRJ   = 498
C      
      IF(CC)WRITE(6,6000)
      IF(MBPT2)WRITE(6,6001)

      CORE_PROJECT   =(Iflags2(h_IFLAGS2_eom_prjct) .GT. 0)
      CORE_SEARCH    =(Iflags2(h_IFLAGS2_ee_search) .EQ. 1)
      DOUBLE_EES     =(Iflags2(h_IFLAGS2_ee_search) .EQ. 4)
C
C  put eigval and oscstr to zero
C
      CALL ZERO(EIGVAL,800)
      CALL ZERO(OSCSTR,800)
      MAXDAV0 = 70
C
C READ IN THE NUMBER OF ROOTS REQUESTED IN EACH SYMMETRY BLOCK
C
      DO 10 IRREPX=1,NIRREP
        ISIDE=1
        IF(NROOT(IRREPX).EQ.0)GOTO 10
        WRITE(6,2000)IRREPX,NROOT(IRREPX)
        CALL NEWLST(IRREPX,SCR,MAXCOR,IUHF)
        IF(NSIZEC.EQ.0)GOTO 10
C
C PUT DIAGONAL PART OF HBAR ON LIST 472
C
        CALL HBARDIAG(IRREPX,SCR,MAXCOR,IUHF)
C
C  DETERMINE EXCITATION PATTERN
C
        DO IOPT = 1,2
          CALL CALCEXCP(IUHF, SCR, MAXCOR, IRREPX,
     &       .TRUE., .FALSE.,NSIZEC,IOPT)
          CALL PUTEXCP(IUHF, SCR, MAXCOR, IRREPX, 
     &       IPROJECT, LISTH0, ICOLPR1, NSIZEC,IOPT)
        ENDDO
C
        N1AA = IRPDPD(IRREPX,ISYTYP(1,19))
        IF (IUHF .NE. 0) THEN
          NDIM = N1AA + IRPDPD(IRREPX,ISYTYP(1,20))
        ELSE
          NDIM = N1AA
        ENDIF
        WRITE(6,2500) NDIM, NSIZEC
C
        NUMROOT = NROOT(IRREPX)
        CALL UPDMOI(NUMROOT,NDIM,IRREPX,LISTRIGHT, 0, 0)
C
C Create pointers for the projection vectors. We need NUMROOT of them
C for each irrep. 

        IF (IPROJECT .GE. 2) CALL UPDMOI(NUMROOT,NDIM,IRREPX,
     &                                   LISTPRJ, 0, 0)
        CALL UPDMOI(NUMROOT,1,IRREPX,LISTE, 0, 0)
        IF (ESPROP) CALL UPDMOI(NUMROOT,NDIM,IRREPX,LISTLEFT, 0, 0)
C 
C At this point we can create projection vectors. They span only
C singles. 
C
        IF (IFLAGS2(120) .EQ. 5) THEN
           Write(6,"(a)") "Excitations energies via NTO-projection."

           CALL DRIVE_NTO_EE_STATE_PRJCT_4PEOM(SCR,MAXCOR,IRREPX,
     &                                         IUHF,NDIM,NUMROOT,
     &                                         LISTPRJ)
          DO ITOP = 1,2
             CALL PUTEXCP(IUHF, ICORE, MAXCOR/IINTFP, IRREPX,
     &                   IPROJECT, LISTH0, ICOLPR1, NSIZEC,ITOP)
          ENDDO

        ENDIF 
C
        CALL GETLST(SCR, 1, NUMROOT, 1, IRREPX, LISTETDA)
C
C  CHECK IF ALL ENERGIES ARE WITHIN A RANGE OF 0.3 A.U. THEN
C  WE DO A BLOCK DAVIDSON, OTHERWHISE WE TURN TO ALTERNATIVE ALGORITHM.        
C
        EMIN = 1000000.0D0
        EMAX = -1000000.0D0
        DO I = 1, NUMROOT
          IF (SCR(I) .LT. EMIN) EMIN = SCR(I)
          IF (SCR(I) .GT. EMAX) EMAX = SCR(I)
        ENDDO
C
C No block Davidson for NTO projection (IFLAGS2(120) .EQ. 5)
C
        DOBLOCK = (EMAX - EMIN) .LE. 0.3
        lockleft = doblock
C
C        write(6,*) ' @doparteom: determine doblock'
C        write(6,*) emin, emax, doblock
C        write(6,*)
C
C Initialize reordering map
C
        DO 195 I = 1, NUMROOT
          IREORD(I) = I
 195    CONTINUE
C
        MAXPOWER = 3
C
C  DETERMINE APPROXIMATE MAXDAVM (THE ABSOLUTE MAXIMUM FOR MAXDAVM)
C
        MXCOR = MAXCOR - 2*NSIZEC - NSIZEC / 4 - 2 * NDIM * NUMROOT
        MAXDAVM = MXCOR / (NDIM * (MAXPOWER + 1))
C
        MAXDAV = MIN(MAXDAV0, NDIM, MAXDAVM)
        IF (TRIPLET) MAXDAV = MIN(MAXDAV, NDIM/2)
C
        IF (DOBLOCK) THEN
C
          CALL TIMER(1)      
C
C DETERMINE A CRUDE ESTIMATE OF EIGENVALUES
C
        WRITE(6,*) ' First determine a crude estimate of eigenvalues'
c
        IZOUT = 1
        IZLEFT = IZOUT + NDIM * NUMROOT
        IEOUT = IZLEFT + NDIM * NUMROOT
        IEOLD = IEOUT + NUMROOT
        IZ = IEOLD + NUMROOT
        IHZ = IZ + NDIM * MAXDAV
        IE0 = IHZ + NDIM * MAXDAV * MAXPOWER
        IH0 = IE0 + MAXDAV
        ISCR = IH0 + NDIM
C
        WRITE(6,*) ' Memory used for Block Davidson ', ISCR
        write(6,*)
        WRITE(6,4000)
        IFIRST = 1
        ILAST = 1
        IF (ESPROP) ILAST = 2
C
        DO 200 ISIDE = IFIRST, ILAST
C
        WRITE(6,*)
        WRITE(6,7000)HAND(ISIDE)
        WRITE(6,*)
C
        DO IROOT = 1, NUMROOT
          ICONV(IROOT) = 1
        ENDDO
C
        IF (ISIDE .EQ. 1) THEN
          CALL GETLST(SCR(IZOUT), 1, NUMROOT, 1, IRREPX, LISTVTDA)
          CALL GETLST(SCR(IEOUT), 1, NUMROOT, 1, IRREPX, LISTETDA)
C
C  SUBTRACT CONSTANT SHIFT FROM TDA EIGENVALUES, AS GIVEN BY EXPERIENCE
C
          TDASHIFT = 0.01D0
          DO I = 1, NUMROOT
            SCR(IEOUT+I-1) = SCR(IEOUT+I-1) - TDASHIFT
          ENDDO
        ELSE
          CALL GETLST(SCR(IZOUT), 1, NUMROOT, 1, IRREPX, LISTLEFT)
          CALL GETLST(SCR(IEOUT), 1, NUMROOT, 1, IRREPX, LISTE)
        ENDIF
C
        TOL = 5.0D-5
C
        SHIFTROOT = .TRUE.
        CALL ED_DAVID(IUHF, IRREPX, ISIDE, NUMROOT, MAXPOWER, MAXDAV,
     &   N1AA, NDIM, SCR(IZOUT), SCR(IZLEFT), SCR(IEOUT),
     &     SCR(IEOLD), SCR(IZ),
     &     SCR(IHZ), SCR(IE0), SCR(IH0), ICONV, SCR(ISCR),
     &     MAXCOR-ISCR+1, 2, TOL, ESPROP .AND. (ISIDE .EQ. 1),
     &     LOCKLEFT) 
C
C  PUT CRUDE ESTIMATES ON LISTRIGHT / LISTLEFT      
C
        IF (ISIDE .EQ. 1) THEN
          CALL PUTLST(SCR(IZOUT), 1, NUMROOT, 1, IRREPX, LISTRIGHT)
          CALL PUTLST(SCR(IEOUT), 1, NUMROOT, 1, IRREPX, LISTE)
          IF (ESPROP) CALL PUTLST(SCR(IZLEFT), 1, NUMROOT, 1,
     &       IRREPX, LISTLEFT)
C
C Sort the energies to get a map from the TDA order to the EOM energy order
C
          CALL PIKSR2(NUMROOT, SCR(IEOUT), IREORD)
        ELSE
          CALL PUTLST(SCR(IZOUT), 1, NUMROOT, 1,IRREPX, LISTLEFT)
        ENDIF
C
  200 CONTINUE
C
      CALL TIMER(1)      
      TBLOCK = TIMENEW
      WRITE(6,*) 
      write(6,1000) TBLOCK
 1000 Format(T3,' Time spent in block Davidson ',F12.4)
      write(6,*)
C
C  PRINT OUT A SUMMARY OF THE CALCULATION UP TILL NOW
C
      IZTDA = 1
      IZRIGHT = IZTDA + NDIM
      IZLEFT = IZRIGHT + NDIM
      WRITE(6,*)
      WRITE(6,2900) IRREPX
      WRITE(6,2800)
      DO IROOT = 1, NUMROOT
        CALL GETLST(SCR(IZTDA),IROOT,1,1,IRREPX,LISTVTDA)
        CALL GETLST(SCR(IZRIGHT),IROOT,1,1,IRREPX,LISTRIGHT)
        IF (ESPROP) CALL GETLST(SCR(IZLEFT),IROOT,1,1,IRREPX,LISTLEFT)
        CALL GETLST(ROOT,IROOT,1,1,IRREPX,LISTE)
        ORIGHT = ABS(SDOT(NDIM, SCR(IZTDA), 1, SCR(IZRIGHT), 1))
        IF (ESPROP) THEN
          OLEFT = ABS(SDOT(NDIM, SCR(IZTDA), 1, SCR(IZLEFT), 1))
        ELSE
          OLEFT = 0.0D0
        ENDIF
        WRITE(6,2950) IROOT, ROOT*FACTEV, ORIGHT, OLEFT
      ENDDO
      WRITE(6,*)  
      IF (NUMROOT .GT. 1) WRITE(6,1030)
C
C  IF NOT DOBLOCK
C
      ELSE
          CALL GETLST(SCR, 1, NUMROOT, 1, IRREPX, LISTVTDA)
          CALL PUTLST(SCR, 1, NUMROOT, 1, IRREPX, LISTRIGHT)
          IF (ESPROP) CALL PUTLST(SCR, 1, NUMROOT, 1, IRREPX,
     &                             LISTLEFT)
          CALL GETLST(SCR, 1, NUMROOT, 1, IRREPX, LISTETDA)
          CALL PUTLST(SCR, 1, NUMROOT, 1, IRREPX, LISTE)
      ENDIF
C
C DETERMINE FULLY CONVERGED STATES
C
        IFIRST = 1
        ILAST = 1
        IF (ESPROP) ILAST = 2
        ITYPECONV = 1
        IF (ESPROP) ITYPECONV = 2
C
C LOOP OVER NUMBER OF DESIRED ROOTS
C SG 2/11/97 First see if user specified which root to follow
C
        CALL CHEXRT(NROOT(IRREPX))
C
        DO 100 IROOT = 1, NROOT(IRREPX)
C
          IRTNEW = IREORD(IROOT)
C
          WRITE(6,1010) IROOT
          WRITE(6,*)
C
          DO 300 ISIDE = IFIRST, ILAST
C
            IF (ISIDE .EQ. 2) THEN
              WRITE(6,*)
              WRITE(6,1020) IROOT
              WRITE(6,*)
            ENDIF
C
            IF (ISIDE .EQ. 1) THEN
              MAXPOWER = 2
            ELSE
              IF (LOCKLEFT) THEN
                MAXPOWER = 1
              ELSE
                MAXPOWER = 2
              ENDIF
            ENDIF
        NUMROOT = 1
        N1AA = IRPDPD(IRREPX,ISYTYP(1,19))
        IF (IUHF .NE. 0) THEN
          NDIM = N1AA + IRPDPD(IRREPX,ISYTYP(1,20))
        ELSE
          NDIM = N1AA
        ENDIF
C
        IZOUT = 1
        IZLEFT = IZOUT + NDIM * NUMROOT
        IEOUT = IZLEFT + NDIM * NUMROOT
        IEOLD = IEOUT + NUMROOT
        IZ = IEOLD + NUMROOT
        IHZ = IZ + NDIM * MAXDAV
        IE0 = IHZ + NDIM * MAXDAV * MAXPOWER
        IH0 = IE0 + MAXDAV
        ISCR = IH0 + NDIM
C
C GET INITIAL GUESS
C
        IF (ISIDE .EQ. 1) THEN
          CALL GETLST(SCR(IZOUT),IRTNEW, NUMROOT,1,IRREPX,LISTRIGHT)
          CALL GETLST(SCR(IEOUT),IRTNEW, NUMROOT,1,IRREPX,LISTE)
        ELSE
          CALL GETLST(SCR(IZOUT),IRTNEW, NUMROOT,1,IRREPX,LISTLEFT)
          SCR(IEOUT) = EIGVAL(IROOT, IRREPX)
        ENDIF
C
        TOL = 10.0D0**(-NTOL)
        SHIFTROOT = .TRUE.
        CALL ED_DAVID(IUHF, IRREPX, ISIDE, NUMROOT, MAXPOWER, MAXDAV,
     &   N1AA, NDIM, SCR(IZOUT), SCR(IZLEFT),SCR(IEOUT),
     &           SCR(IEOLD), SCR(IZ),
     &     SCR(IHZ), SCR(IE0), SCR(IH0), ICONV, SCR(ISCR),
     &     MAXCOR-ISCR+1, ITYPECONV, TOL, .FALSE., LOCKLEFT) 
C
C Z and root are determined. at this point other properties can
C be calculated if desired. 
C
            ROOT = SCR(IEOUT)
            IF (ISIDE .EQ. 1) EIGVAL(IROOT, IRREPX) = ROOT
C
            CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ECC)
            IF(CC)THEN
              WRITE(6,99)ECC+ROOT
            ELSE
              WRITE(6,199)ECC+ROOT
            ENDIF
   99       FORMAT(T3,' Total P-EOM-CCSD electronic energy ',
     &         F20.12,' a.u.')
  199       FORMAT(T3,' Total P-EOM-MBPT(2) electronic energy ',
     &         F20.12,' a.u.')
            CALL PUTREC(20,'JOBARC','TOTENER2',IINTFP,ECC+ROOT)
C
C CALCULATE TOTAL CONVERGED VECTOR
C
            CALL PUTLST(SCR(IZOUT), 1, 1, 1, 1, 490)
            IF (IUHF .NE. 0) CALL PUTLST(SCR(IZOUT+N1AA),1,1,1,2,490)
C
            I000 = 1
            I010 = I000 + NSIZEC
            I020 = I010 + NSIZEC
C
C CALCULATE DOUBLES PART
C
            SHIFTROOT = .FALSE.
            CALL PHBARXC(SCR, MAXCOR*IINTFP,IUHF, ISIDE, IRREPX,
     &         .FALSE., SD, .FALSE., .FALSE., ROOT, 1)       
C
            CALL LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,490,0,443,NSIZEC,
     &         .FALSE.)
C
C  print out main components of eigenvector
C
            CALL SCOPY(NSIZEC,SCR(I000),1,SCR(I010),1)
            IF(IUHF.EQ.0)THEN
              CALL SPNTSING(IRREPX,SCR(I010),SCR(I020),MAXCOR-I010+1)
            ENDIF
            Z=SDOT(NSIZEC,SCR(I000),1,SCR(I010),1)
            X=ONE/SQRT(Z)
            CALL SSCAL (NSIZEC,X,SCR(I000),1)
C
C PRINT MAIN COMPONENTS OF EIGENVECTOR
C
            CALL UPDATES(IRREPX,SCR(I000),444,0,490,IUHF)
            CALL PRMAINX(IUHF, SCR(I010), MAXCOR, IRREPX,IROOT)
C
            IF(ESPROP)THEN
              IF(ISIDE.EQ.1)THEN
                CALL RNORM(IRREPX,NSIZEC,ROOT,SCR(I000),SCR(I010),
     &             MAXCOR-I010+1,IUHF,Z0)
                CALL PUTLST(SCR(I000),1+ISIDE,1,1,1,472)
                IF(LAMTHERE)THEN
                  CALL TDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &               (ONE+DFLOAT(1-IUHF)),Z0,ROOT,CORE_SEARCH,
     &               DOUBLE_EES)
                ENDIF
              ELSEIF(ISIDE.EQ.2)THEN
                IF (ABS(ROOT - EIGVAL(IROOT,IRREPX)) .GT. 1.E-3) THEN
                  write(6,*) ' left and right eigenvectors do not match'
                  OSCSTR(IROOT,IRREPX) = 1.0D20
                ELSE
                  CALL PUTLST(SCR(I000),1+ISIDE,1,1,1,472)
                  CALL UPDATES(IRREPX,SCR(I000),444,0,490,IUHF)
                  CALL LNORM(IRREPX,ROOT,SCR,MAXCOR,IUHF,Z,Z0)
                  IF(LAMTHERE)THEN
                    CALL TDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &                 Z*(ONE+DFLOAT(1-IUHF)),Z0,ROOT,CORE_SEARCH,
     &                 DOUBLE_EES)
                    CALL PRINTSUM(ROOT,FSTR,CORE_SEARCH,DOUBLE_EES)
                    OSCSTR(IROOT,IRREPX) = FSTR
                  ENDIF
                  CALL GETLST(SCR(I000),3,1,1,1,472)
                  CALL SSCAL (NSIZEC,Z,SCR(I000),1)
                  CALL PUTLST(SCR(I000),3,1,1,1,472)
                  IF (IFLAGS(91) .LT. 2) THEN
                    CALL EDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &                 Z*(ONE+DFLOAT(1-IUHF)),Z0)
                  ELSEIF (IROOT .EQ. NROOT(IRREPX)) THEN
                    call inigam(iuhf)
                    CALL EDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &                 Z*(ONE+DFLOAT(1-IUHF)),Z0)
C
C CALCULATE EOM-CCSD TWO-PARTICLE DENSITY
C
                    CALL TPDENS(SCR,IINTFP*MAXCOR,IUHF,Z0)
                    call aces_fin
                    STOP
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
C
  300     CONTINUE
C
  100   CONTINUE
C
   10 CONTINUE
C
      TOTROOT = 0
      DO I = 1, NIRREP
        TOTROOT = TOTROOT + NROOT(I)
      ENDDO
      IF (TOTROOT .GT. 1) CALL TABLEE(EIGVAL,EIGVAL_T,OSCSTR,NATURE,
     &                                BGN,BGN_IRP,END,END_IRP,
     &                                IFLAGS(87))
C
      RETURN
C
 1010 FORMAT(T3,'Solving for root # ', I3, '.')
 1020 FORMAT(T3,'Solving for left hand root # ', I3, '.')
 1030 FORMAT(T3,'@DOPARTEOM-I, Excited states will be calculated in',
     &   ' EOM energy order.')
 2000 FORMAT(/,T3,'Beginning symmetry block ',I3,'.',I4,
     &   ' roots requested.')
 2500 FORMAT(T3, 'Dimension singles ',I4, ' Total Dimension ',I8)
 2800 FORMAT(T3,'root #',T12,' eigenvalue (eV)', T33,
     &   ' overlap TDA solution ',/, T36, ' Right', T47, 'Left')
 2900 FORMAT(T3, 'Summary of Block-Davidson Results for irrep',I3)
 2950 FORMAT(T3,I4,T10,F15.8,T33,F8.3,T43,F8.3)
 4000 FORMAT(T18,'Subspace',T36,'Eigenvalue',T51,' ',/,
     &       T8,'ROOT #',
     &       T18,'Dimension',T31,'(a.u.)',T46,'(eV)',T59,
     &       'Residual')
 6000 FORMAT(T3,'@DOPARTEOM-I, Excitation energies computed by the ',
     &          'P-EOM-CCSD method.')
 6001 FORMAT(T3,'@DOPARTEOM-I, Excitation energies computed by the ',
     &          'P-EOM-MBPT(2) method.')
 7000 FORMAT(T3,'Block_Davidson:, ',A,'-hand eigenvectors will ',
     &   'be computed.')
C
      END
