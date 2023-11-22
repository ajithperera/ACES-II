      SUBROUTINE GET_CORR_DENS(DENS, SCR1, SCR2, SCRA, SCRB, SCRTMP,
     &                         SCRN, SCR3, SCR, LDIM1, LDIM2, NBAS, 
     &                         NBASX, IUHF)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
#include "machsp.com"
#include "symm2.com"
C
      DIMENSION DENS((IUHF+1)*LDIM1), SCR1(NBAS,NBAS), SCR2(LDIM2),
     &          SCR3(NBAS,NBAS), SCRA(LDIM2), SCRB(LDIM2),
     &          SCRN(NBAS,NBAS), SCRTMP(LDIM2), SCR(NBAS,NBAS)
      DATA ZILCH, ONE / 0.0D0, 1.0D0/
C
      ZTEST = 0.0D+00
      CALL ZERO(SCR1,NBAS*NBAS)
      CALL ZERO(SCRN,NBAS*NBAS)
C
      CALL GETREC(20,'JOBARC','RELDENSA',NBAS*NBAS*IINTFP,SCR1)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The Alpha Density read"
      call output(SCR1, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
C
      CALL EIG(SCR1,SCRN,1,NBAS,-1)
C
#ifdef _REORDER
      ZTEST=SNRM2(NBAS*NBAS,SCRN,1)
      IF(ZTEST.NE.ZILCH) THEN
         CALL REOREV(SCRN,SCR3,NBAS,1)
C
#ifdef _DEBUG_LVL0
        Write(6,"(a)") "Reordered orbs: Alpha"
        call output(SCRN, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
       DO JBF = 1, NBAS
           Z=SNRM2(NBAS, SCRN(1,JBF), 1)
           CALL SSCAL(NBAS,ONE/Z, SCRN(1,JBF),1)
       ENDDO
       CALL SORTHO(SCRN,SCR3,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
      ENDIF
#endif
C
      CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR3)
      CALL XGEMM('N', 'N', NBAS, NBAS, NBAS, ONE, SCR3, NBAS,
     &           SCRN, NBAS, ZILCH, SCR, NBAS)

#ifdef _DEBUG_LVLM1
      Write(6, "(a)") "The Alpha natural orbitals"
      call output(buf(ioffb2), 1, nbas, 1, nbas, nbas, nbas, 1)
#endif
      CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR)

#ifdef _DEBUG_LVL0
      Write(6,*) "The occupation number"
      call output(SCR1, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
      CALL MKDDEN_FULL(SCRN,SCR2,SCR3,SCR1,NBAS,IUHF)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The Alpha Density after MKDDEN"
      call output(SCR3, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
C
      CALL MO2AO2(SCR3, SCRN, SCR1, SCR, NBAS, 1)
C
      DO I=1,NIRREP
          CALL GETBLK(SCRN,SCR2,NBFIRR(I),NBAS,IREPS(I))
          CALL SQUEZ2(SCR2,DENS(ITRIOF(I)),NBFIRR(I))
      ENDDO
C
      IF (IUHF.EQ.1) THEN
         ZTEST = 0.0D+00
         CALL ZERO(SCR1,NBAS*NBAS)
         CALL ZERO(SCRN,NBAS*NBAS)
         CALL GETREC(20,'JOBARC','RELDENSB',NBAS*NBAS*IINTFP,SCR1)
#ifdef _DEBUG_LVLM1
         Write(6,*) "The Beta Density read"
         call output(SCR1, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
         CALL EIG(SCR1,SCRN,1,NBAS,-1)
C
#ifdef _REORDER
         ZTEST=SNRM2(NBAS*NBAS,SCRN,1)
         IF (ZTEST.NE.ZILCH) THEN
            CALL REOREV(SCRN,SCR3,NBAS,2)

            DO JBF = 1, NBAS
               Z=SNRM2(NBAS, SCRN(1,JBF), 1)
               CALL SSCAL(NBAS,ONE/Z, SCRN(1,JBF),1)
            ENDDO
         CALL SORTHO(SCRN,SCR3,SCRA,SCRB,SCRTMP,LDIM2,NBAS)
      ENDIF
#endif
C
      CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,SCR3)
      CALL XGEMM('N', 'N', NBAS, NBAS, NBAS, ONE, SCR3, NBAS,
     &           SCRN, NBAS, ZILCH, SCR, NBAS)

#ifdef _DEBUG_LVLM1
      Write(6, "(a)") "The Alpha natural orbitals"
      call output(buf(ioffb2), 1, nbas, 1, nbas, nbas, nbas, 1)
#endif
      CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR)


         CALL MKDDEN_FULL(SCRN,SCR2,SCR3,SCR1,NBAS,IUHF)
C
#ifdef _DEBUG_LVLM1
         Write(6,*) "The Beta Density after MKDDEN"
         call output(SCR3, 1, NBAS, 1, NBAS, NBAS, NBAS, 1)
#endif
         CALL MO2AO2(SCR3, SCRN, SCR1, SCR, NBAS, 2)
C
         DO I=1,NIRREP
            CALL GETBLK(SCRN,SCR2,NBFIRR(I),NBAS,IREPS(I))
            CALL SQUEZ2(SCR2,DENS(LDIM1+ITRIOF(I)),NBFIRR(I))
         ENDDO  
C
      ENDIF

      RETURN
      END
