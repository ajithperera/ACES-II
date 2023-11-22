#include "flags.h"
      SUBROUTINE DOBWPT2(IUHF, SCR, MAXCOR)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION ICONV(100)
      DIMENSION SCR(MAXCOR)
      LOGICAL ESPROP
      LOGICAL LAMTHERE, LOCKLEFT, SHIFTROOT
      LOGICAL MBPT2, CC, SS, SD, DS, DD, CCD,RCCD,DRCCD
      LOGICAL LCCD,LCCSD
      LOGICAL CORE_PROJECT,CORE_SEARCH,DOUBLE_EES
      INTEGER DIRPRD
      CHARACTER*5 HAND(2)
C
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREP0(255,2),DIRPRD(8,8)
      COMMON/CALCINFO/NROOT(8)
      COMMON/ROOTS/EIGVAL(100,8), OSCSTR(100,8)
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
      COMMON/TIMSUB/TDAVID, TMULT
      COMMON/SHFTROOT/SHIFTROOT
C
      DATA ONE /1.0D0/
      DATA HAND /'right',' left'/
C
      LOCKLEFT = (IPATTERN .NE. 1)
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
        MAXPOWER = 3
        MAXDAV = 1
        NUMROOT = 1
        IFIRST = 1
        ILAST = 1
        IF (ESPROP) ILAST = 2
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
        WRITE(6,4000)
        ISIDE = 1
        WRITE(6,*)
        WRITE(6,7000)HAND(ISIDE)
        WRITE(6,*)
C
        DO 200 IROOT = 1, NROOT(IRREPX)
C
          ICONV(1) = 1
          ISIDE = 1
C
        IF (ISIDE .EQ. 1) THEN
          CALL GETLST(SCR(IZOUT), IROOT, NUMROOT, 1, IRREPX, LISTVTDA)
          CALL GETLST(SCR(IEOUT), IROOT, NUMROOT, 1, IRREPX, LISTETDA)
C
C  SUBTRACT CONSTANT SHIFT FROM TDA EIGENVALUES, AS GIVEN BY EXPERIENCE
C
          TDASHIFT = 0.01D0
          DO I = 1, NUMROOT
            SCR(IEOUT+I-1) = SCR(IEOUT+I-1) - TDASHIFT
          ENDDO
        ENDIF
C
        TOL = 5.0D-5
C
        SHIFTROOT = .TRUE.
        CALL ED_DAVID(IUHF, IRREPX, ISIDE, NUMROOT, MAXPOWER, MAXDAV,
     &   N1AA, NDIM, SCR(IZOUT), SCR(IZLEFT), SCR(IEOUT),
     &     SCR(IEOLD), SCR(IZ),
     &     SCR(IHZ), SCR(IE0), SCR(IH0), ICONV, SCR(ISCR),
     &     MAXCOR-ISCR+1, 2, TOL, .FALSE.,
     &     LOCKLEFT) 
C
C Z and root are determined. at this point other properties can
C be calculated if desired. 
C
        CALL PUTLST(SCR(IEOUT), IROOT, NUMROOT, 1, IRREPX, LISTETDA)
            ROOT = SCR(IEOUT)
            IF (ISIDE .EQ. 1) EIGVAL(IROOT, IRREPX) = ROOT
C
            IF (IFLAGS(87) .EQ. 8) THEN
              CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ECC)
              IF(CC)THEN
                WRITE(6,99)ECC+ROOT
              ELSE
                WRITE(6,199)ECC+ROOT
              ENDIF
   99         FORMAT(T3,' Total EOM-BW-CCSD electronic energy ',
     &           F20.12,' a.u.')
  199         FORMAT(T3,' Total EOM-BWPT(2) electronic energy ',
     &           F20.12,' a.u.')
              CALL PUTREC(20,'JOBARC','TOTENER2',IINTFP,ECC+ROOT)
            ENDIF
C
            DO 50 ISIDE = IFIRST, ILAST
C
C CALCULATE TOTAL CONVERGED VECTOR
C
            CALL GETLST(SCR(IZOUT), IROOT, NUMROOT, 1, IRREPX, LISTVTDA)
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
            IF (IFLAGS(87) .EQ. 8) THEN
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
              IF (ISIDE .EQ. 1) THEN
                WRITE(6,*) '   Right hand eigenvector '
              ELSE
                WRITE(6,*) '   Left hand eigenvector '
              ENDIF
              CALL UPDATES(IRREPX,SCR(I000),444,0,490,IUHF)
              CALL PRMAINX(IUHF, SCR(I010), MAXCOR, IRREPX,IROOT)
            ENDIF
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
                    IF (IFLAGS(87) .EQ. 8) THEN
                      CALL PRINTSUM(ROOT,FSTR,CORE_SEARCH,DOUBLE_EES)
                    ELSE
                      CALL CALCFSTR(ROOT,FSTR)
                    ENDIF
                    OSCSTR(IROOT,IRREPX) = FSTR
                  ENDIF
                  CALL GETLST(SCR(I000),3,1,1,1,472)
                  CALL SSCAL (NSIZEC,Z,SCR(I000),1)
                  CALL PUTLST(SCR(I000),3,1,1,1,472)
                  IF (IFLAGS(87) .EQ. 8) THEN
                    IF (IFLAGS(91) .LT. 2) THEN
                      CALL EDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &                   Z*(ONE+DFLOAT(1-IUHF)),Z0)
                    ELSEIF (IROOT .EQ. NROOT(IRREPX)) THEN
                      call inigam(iuhf)
                      CALL EDENS(IRREPX,SCR,MAXCOR,IUHF,ISIDE,
     &                   Z*(ONE+DFLOAT(1-IUHF)),Z0)
C
C CALCULATE EOM-CCSD TWO-PARTICLE DENSITY
C SG 3/96 two-particle density is not available for eom-bwpt2
                      CALL TPDENS(SCR,IINTFP*MAXCOR,IUHF,Z0)
                      call aces_fin
                      STOP
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
C
   50     CONTINUE
C
  200     CONTINUE
C
  10   CONTINUE
C
      TOTROOT = 0
      DO I = 1, NIRREP
        TOTROOT = TOTROOT + NROOT(I)
      ENDDO
C
      IF (TOTROOT .GT. 1) CALL TABLEE(EIGVAL,OSCSTR,8)
C
      WRITE(6,*)
      write(6,1000) TMULT
 1000 Format(T3,' Time spent in block Davidson ',F12.4)
      WRITE(6,*)
      RETURN
C
 1010 FORMAT(T3,'Solving for root # ', I3, '.')
 1020 FORMAT(T3,'Solving for left hand root # ', I3, '.')
2000  FORMAT(/,T3,'Beginning symmetry block ',I3,'.',I4,
     &   ' roots requested.')
 2500 FORMAT(T3, 'Dimension singles ',I4, ' Total Dimension ',I8)
6000  FORMAT(T3,'@DOBWPT2-I, Excitation energies from the ',
     &          'Brillouin-Wigner EOM-CCSD method.')
6001  FORMAT(T3,'@DOBWPT2-I, Excitation energies by the ',
     &   'Brillouin-Wigner EOM-MBPT(2) method.')
7000  FORMAT(T3,'Block_Davidson:, ',A,'-hand eigenvectors will ',
     &   'be computed.')
4000  FORMAT(T18,'Subspace',T36,'Eigenvalue',T51,' ',/,
     &       T8,'ROOT #',
     &       T18,'Dimension',T31,'(a.u.)',T46,'(eV)',T59,
     &   'Residual')
 2900 FORMAT(T3, 'Summary of Block-Davidson Results for irrep',I3)
 2800 FORMAT(T3,'root #',T12,' eigenvalue (eV)', T33,
     &   ' overlap TDA solution ',/, T36, ' Right', T47, 'Left')
 2950  FORMAT(T3,I4,T10,F15.8,T33,F8.3,T43,F8.3)
C
      END
