
      SUBROUTINE EXPORT_GUESS_MO_ORBS(SCR1, SCRN, SCRA, SCRB,
     &                                SCR2, LDIM2, NBAS, IUHF, ROHF)

       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
       LOGICAL GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR,
     &         ROHFMO,FOLLOW,ROHF
       INTEGER SWAP,READMO,WRITMO
C
#include "machsp.com"
#include "symm2.com"
C
       DIMENSION SCR1(NBAS*NBAS),SCRN(NBAS*NBAS),
     &           SCRA(LDIM2),SCRB(LDIM2),SCR2(LDIM2)

      COMMON /POPUL/  NOCC(16)
      COMMON /GSCOMA/ GSSOPT,GSSALW,GSSALT,GSSLOK,GSSRED,GSSWRT,GSSUFR
      COMMON /GSCOMB/ SWAP(4,8,2),LOCK(8,2),IPRTGS(8,2),ISTOP1,ISTOP2,
     &                READMO,WRITMO,IUHFRHF,LUGSS
      COMMON /GSCOMC/ ROHFMO
      COMMON /FILES/  LUOUT,MOINTS
      COMMON /FLAGS/  IFLAGS(100)
C
C     Initialize the alpha and beta full MO coefficient matrices.
C
       CALL ZERO(SCR1,NBAS*NBAS)
       CALL ZERO(SCRN,NBAS*NBAS)
C
       DO  450 I=1,NIRREP
       IF(NBFIRR(I).EQ.0) GOTO 450
C
       IF(GSSRED)THEN
C
        WRITE(6,6030)
 6030   FORMAT(' @INITGES-I, Initial guess orbitals read from OLDMOS. ')
C
C     Read orbitals from OLDMOS for this symmetry block.
C
        OPEN(LUGSS,FILE='OLDMOS',STATUS='OLD',ACCESS='SEQUENTIAL',
     1       FORM='FORMATTED')
        CALL GETMOS(SCRA,NBFIRR(I),1,I,NIRREP,NBFIRR,LUGSS)
C
        IF(IUHF.GT.0.AND..NOT.ROHF)THEN
         REWIND LUGSS
         IF(GSSUFR)THEN
C     Read first set of MOs on OLDMOS into beta array.
          CALL GETMOS(SCRB,NBFIRR(I),1,I,NIRREP,NBFIRR,LUGSS)
         ELSE
C     Read second set of MOs on OLDMOS into beta array.
          CALL GETMOS(SCRB,NBFIRR(I),2,I,NIRREP,NBFIRR,LUGSS)
         ENDIF
        ELSE
c YAU : old
c        CALL ICOPY(NBFIRR(I)*NBFIRR(I)*IINTFP,SCRA,1,SCRB,1)
c YAU : new
         CALL DCOPY(NBFIRR(I)*NBFIRR(I),SCRA,1,SCRB,1)
c YAU : end
        ENDIF
C
        CLOSE(LUGSS,STATUS='KEEP')
C
C
       ENDIF
C
C     The initial guess orbitals have now been generated. We now deal
C     with issues such as printing and swapping.
C
C     If this has been requested, print the initial guess.
C
       IF(IPRTGS(I,1).GT.0)THEN
        WRITE(LUOUT,2010) I
        CALL OUTPUT(SCRA,1,NBFIRR(I),1,NBFIRR(I),NBFIRR(I),NBFIRR(I),1)
       ENDIF
C
       IF(IPRTGS(I,2).GT.0)THEN
        WRITE(LUOUT,2020) I
        CALL OUTPUT(SCRB,1,NBFIRR(I),1,NBFIRR(I),NBFIRR(I),NBFIRR(I),1)
       ENDIF
 2010  FORMAT(' @INITGES-I, Alpha initial guess. Symmetry block ',I3)
 2020  FORMAT(' @INITGES-I, Beta  initial guess. Symmetry block ',I3)
C
C     Construct alpha density matrix.
C
C     First swap orbitals (if this has been requested)
C
       IF(GSSALT)THEN
        DO  410 IPAIR=1,2
        IF(SWAP(2*IPAIR-1,I,1).GT.0)THEN
         I1 = SWAP(2*IPAIR-1,I,1)
         I2 = SWAP(2*IPAIR  ,I,1)
         WRITE(6,2030) I,I1,I2
 2030    FORMAT(' @INITGES-I, Alpha spin, symmetry block ',I3,' . ',/,
     &   '             Orbitals ',I4,' and ',I4,' will be swapped. ')
         CALL SWAPORB(SCRA,SCR2,NBFIRR(I),I1,I2)
        ENDIF
  410   CONTINUE
       ENDIF
C
C     Form the contribution of this symmetry block of orbitals to the
C     whole coefficient matrix. (Alpha).
C
       CALL PUTBLK(SCR1,SCRA,NBFIRR(I),NBAS,IREPS(I))
C
C     Construct beta density matrix if this is not closed-shell RHF.
cjp bug ... forgotten iuhf test, scr1 gets incorrectly overwritten
cjp
      if(iuhf.eq.1) then
cjp
       IF(GSSALT)THEN
        DO  420 IPAIR=1,2
        IF(SWAP(2*IPAIR-1,I,2).GT.0)THEN
         I1 = SWAP(2*IPAIR-1,I,2)
         I2 = SWAP(2*IPAIR  ,I,2)
         WRITE(6,2040) I,I1,I2
 2040    FORMAT(' @INITGES-I, Beta  spin, symmetry block ',I3,' . ',/,
     &    '             Orbitals ',I4,' and ',I4,' will be swapped. ')
         CALL SWAPORB(SCRB,SCR2,NBFIRR(I),I1,I2)
        ENDIF
  420   CONTINUE
       ENDIF
C
       CALL PUTBLK(SCRN,SCRB,NBFIRR(I),NBAS,IREPS(I))
cjp
      endif
cjp
C
  450  CONTINUE
C
C     If all we wanted to do was print the guess, stop here.
C
       IF(ISTOP1.GT.0)THEN
        CALL ERREX
       ENDIF
C
C     Dump initial eigenvectors (SO basis) to disk.
C
CSSs       CALL REOREV(SCR1,SCRN, NBAS, 1)
       ZTEST = DNRM2(NBAS*NBAS, SCR1, 1)
       CALL SORTHO(SCR1, SCRN, SCRA, SCRB, SCR2, LDIM2, NBAS)
       ZTEST = DNRM2(NBAS*NBAS, SCR1, 1)
       CALL SORTHO(SCR1, SCRN, SCRA, SCRB, SCR2, LDIM2, NBAS)
       Write(6,*) "REORDERD NATURAL ORBITALS", ZTEST
       call output(SCR1, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
       CALL PUTREC(-1,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
       IF(IUHF.GT.0)THEN
CSSS        CALL REOREV(SCRN,SCR1, NBAS, 2)
        CALL PUTREC(-1,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCRN)
       ENDIF
C
      RETURN
      END

