
      SUBROUTINE EIGVECTORS_FRM_JOBARC(SCR1, SCRN, SCRA, SCRB,
     &                                 SCRTMP, LDIM2, NBAS,
     &                                 IUHF)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
#include "machsp.com"
#include "symm2.com"
C
      DIMENSION SCR1(NBAS*NBAS), SCRN(NBAS*NBAS), SCRA(LDIM2),
     &          SCRB(LDIM2),SCRTMP(LDIM2)
C
      COMMON /FLAGS/IFLAGS(100)
      COMMON /FLAGS2/IFLAGS2(500)
      DATA ZILCH /0.0/
C
C Not Bruckner and HF Stability=on
C
CSSS      IF(IFLAGS(22) .EQ. 0 .AND. IFLAGS(74) .LT. 2)THEN
C I am changing this to IFLAGS(74) .EQ. 1. The SCFEVCAS and
C SCFEVCA0 have same data (within the confines of SCF). The 
C decision has to be made  whether we are going to read 
C SCFEVCA0 or SCFEVCAS. The diffrences is that it is possible
C that SCFEVCA0 can be reorderd (Brueckner?). The
C HFSATBILTY {ON, FOLLOW} reorder the vectors and so 
C read the SCFEVCA0. Also, allow the external program to
C write SCFEVCA0 because that is the record that people are
C familiar with (Localization programs). My feeling is that
C it is better we always read SCFEVCA0. Ajith Perera, 04/2014.
C The IFLAGS(178) correspond to OO-METHOD flag. 

      IF(IFLAGS(22) .EQ. 0 .AND. IFLAGS(74) .LT. 2 .AND.
     &                         IFLAGS2(178) .LE. 0) THEN
C
        ZTEST = 0.0D+00
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCRN,NBAS*NBAS)
        CALL GETREC(-1,'JOBARC','SCFEVCAS',NBAS*NBAS*IINTFP,SCRN)
#ifdef _DEBUG_LVL0
        Write(6,"(2x,a)") "The Alpha eigenvectors read from JOBARC"
        call checksum("SCF_VECA", SCRN, NBAS*NBAS)
CSSS        call output(Scrn, 1, Nbas, 1, Nbas, Nbas, Nbas, 1)
#endif
        ZTEST=SNRM2(NBAS*NBAS,SCRN,1)
        IF(ZTEST.NE.ZILCH) THEN
         CALL SORTHO(SCRN,SCR1,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
        ENDIF
        CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCRN)
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCRN,NBAS*NBAS)
C
        IF(IUHF.EQ.1) THEN
         ZTEST = 0.0D+00
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
         CALL GETREC(-1,'JOBARC','SCFEVCBS',NBAS*NBAS*IINTFP,SCRN)
         ZTEST=SNRM2(NBAS*NBAS,SCRN,1)
         IF(ZTEST.NE.ZILCH) THEN
          CALL SORTHO(SCRN,SCR1,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
         ENDIF
         CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCRN)
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
        ENDIF
C
       ELSE

        ZTEST = 0.0D+00
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCRN,NBAS*NBAS)
        CALL GETREC(-1,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCRN)
        ZTEST=SNRM2(NBAS*NBAS,SCRN,1)

#ifdef _DEBUG_LVL0
        Write(6,"(2x,a)") "The Alpha eigenvectors read from JOBARC"
        call checksum("SCF_VECA0", SCRN, NBAS*NBAS)
C        call output(Scrn, 1, Nbas, 1, Nbas, Nbas, Nbas, 1)
C        Write(*,"(A,F10.6)") "NORM-CHECK =", ZTEST
#endif
        IF(ZTEST.NE.ZILCH) THEN
           CALL REOREV(SCRN,SCR1,NBAS,1)
           CALL SORTHO(SCR1,SCRN,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
        ENDIF
        CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
        CALL ZERO(SCR1,NBAS*NBAS)
        CALL ZERO(SCRN,NBAS*NBAS)
C
        IF(IUHF.EQ.1) THEN
         ZTEST = 0.0D+00
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
         CALL GETREC(-1,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCRN)
         ZTEST=SNRM2(NBAS*NBAS,SCRN,1)
#ifdef _DEBUG_LVL0
        Write(6,"(2x,a)") "The Beta eigenvectors read from JOBARC"
        call checksum("SCF_VECB", SCRN, NBAS*NBAS)
C        Write(*,"(A,F10.6)") "NORM-CHECK =", ZTEST
C        call output(Scrn, 1, Nbas, 1, Nbas, Nbas, Nbas, 1)
#endif 
         IF(ZTEST.NE.ZILCH) THEN
         CALL REOREV(SCRN,SCR1,NBAS,2)
         CALL SORTHO(SCR1,SCRN,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
         ENDIF
         CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCR1)
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
        ENDIF
      ENDIF

      RETURN
      END

