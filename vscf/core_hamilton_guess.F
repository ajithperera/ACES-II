
      SUBROUTINE CORE_HAMILTON_GUESS(SCR1, SCRN, ONEH, SCRA,
     &                                SCR2, SCRB, XFORM, SCRTMP,
     &                                LDIM1, LDIM2, IUHF, NBAS)

       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
#include "molcas.com"
#include "symm2.com"
#include "machsp.com"
C
       DIMENSION SCR1(NBAS*NBAS),SCRN(NBAS*NBAS),ONEH(LDIM1),
     &           XFORM(LDIM1),SCR2(LDIM2),SCRTMP(LDIM2),
     &           SCRA(LDIM2),SCRB(LDIM2)
C
       CALL ZERO(SCR1,NBAS*NBAS)
       CALL ZERO(SCRN,NBAS*NBAS)
C
       DO 110 I=1,NIRREP
       IF(NBFIRR(I).EQ.0) GOTO 110
C
       CALL EXPND2(ONEH(ITRIOF(I)),SCR2,NBFIRR(I))
#ifdef _DEBUG_LVLM
       call output(scr2,1,NBFIRR(I),1,NBFIRR(I),NBFIRR(I),
     +             NBFIRR(I),1)
#endif
       CALL EXPND2(XFORM(ITRIOF(I)),SCRA,NBFIRR(I))
       CALL TRANSP(SCRA,SCRB,NBFIRR(I),NBFIRR(I))
       CALL MXM(SCRB,NBFIRR(I),SCR2,NBFIRR(I),SCRTMP,NBFIRR(I))
       CALL MXM(SCRTMP,NBFIRR(I),SCRA,NBFIRR(I),SCR2,NBFIRR(I))
       CALL ZERO(SCRB,NBFIRR(I)*NBFIRR(I))
       CALL EIG(SCR2,SCRB,NBFIRR(I),NBFIRR(I),0)
C
C     SCRB has core Hamiltonian eigenvectors in the "canonical" basis.
C     Transform these to the SO basis.
C
       CALL MXM(SCRA,NBFIRR(I),SCRB,NBFIRR(I),SCR2,NBFIRR(I))
C
C     Core Hamiltonian orbitals in the SO basis for this symmetry block
C     are in SCR2. Copy into SCRA and SCRB (alpha and beta sets).
C
c YAU : old
c      CALL ICOPY(NBFIRR(I)*NBFIRR(I)*IINTFP,SCR2,1,SCRA,1)
c      CALL ICOPY(NBFIRR(I)*NBFIRR(I)*IINTFP,SCR2,1,SCRB,1)
c YAU : new
       CALL DCOPY(NBFIRR(I)*NBFIRR(I),SCR2,1,SCRA,1)
       CALL DCOPY(NBFIRR(I)*NBFIRR(I),SCR2,1,SCRB,1)
c YAU : end
C
       CALL PUTBLK(SCR1,SCRA,NBFIRR(I),NBAS,IREPS(I))
       IF(IUHF.GT.0)THEN

        CALL PUTBLK(SCRN,SCRB,NBFIRR(I),NBAS,IREPS(I))
       ENDIF
C
  110  CONTINUE
C
#ifdef _DEBUG_LVLM
       Call checksum("GUESS CA",SCR1, NBAS*NBAS)
       Call checksum("GUESS CB",SCRN, NBAS*NBAS)
#endif
       CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
       CALL PUTREC(20,'JOBARC','GUESSVA0',NBAS*NBAS*IINTFP,SCR1)
       IF(IUHF.GT.0)THEN
        CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCRN)
        CALL PUTREC(20,'JOBARC','GUESSVB0',NBAS*NBAS*IINTFP,SCRN)
       ENDIF
C
      RETURN
      END

