
      SUBROUTINE MAKE_GUESS_DEN(SCR1, SCR2, SCRA, SCRN, DENS, NBAS,
     &                          LDIM1, LDIM2, IUHF, ROHF, ROHFMO,
     &                          GUESS_TYPE)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
#include "symm2.com"
#include "machsp.com"
#include "maxbasfn.par"
C
      CHARACTER*11 GUESS_TYPE
      LOGICAL ROHF, ROHFMO
C
      DIMENSION SCR1(NBAS*NBAS),SCRN(NBAS*NBAS), SCR2(LDIM2),
     &          SCRA(LDIM2), DENS((IUHF+1)*LDIM1), DOCC(MAXBASFN*2),
     &          NSUM(16)
C
      COMMON /FILES/  LUOUT,MOINTS
      COMMON /FLAGS/  IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /POPUL/  NOCC(16)
C
      IF (GUESS_TYPE.NE."NDDO"        .AND.
     &    GUESS_TYPE.NE."PROJ_STO-3G" .AND.
     &    GUESS_TYPE.NE."EHT"         .AND. 
     &    GUESS_TYPE.NE. "EXPORT_DENS") THEN
C
           CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
#ifdef _DEBUG_LVLM
          Write(6,*) "@-MAKE_GUESS_DEN: The MOs"
          call  output(scr1, 1, nbas, 1, nbas, nbas, nbas, 1)
#endif 
#ifdef _DPOCC
           CALL MKDP_OCCNOS(DOCC, NOCC, NSUM, NBAS, IUHF, .TRUE.)
#endif
           IROFF = 0
           DO 810 I=1,NIRREP
              CALL GETBLK(SCR1,SCR2,NBFIRR(I),NBAS,IREPS(I))
#ifdef _DPOCC
              IOFF = IROFF + 1 
              CALL MKDDEN_4IRREP(SCR2,SCRA,SCRN,NSUM(I),DOCC(IOFF),
     &                           NBFIRR(I),1,IUHF)
              IROFF = IROFF + NBFIRR(I)
#else
              CALL MKDEN(SCR2,SCRA,NOCC(I),NBFIRR(I),IUHF)
#endif
              CALL SQUEZ2(SCRA,DENS(ITRIOF(I)),NBFIRR(I))
  810      CONTINUE
C
           IF (IUHF .EQ. 1) THEN
               CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,
     &                     SCR1)
               IROFF = 0
               DO 820 I=1,NIRREP
                  IOFF = IROFF  + NBAS + 1
                  CALL GETBLK(SCR1,SCR2,NBFIRR(I),NBAS,IREPS(I))
#ifdef _DPOCC
                  CALL MKDDEN_4IRREP(SCR2,SCRA,SCRN,NSUM(8+I),
     &                               DOCC(IOFF),NBFIRR(I),2,IUHF)
                  IROFF = IROFF + NBFIRR(I)
#else
                  CALL MKDEN(SCR2,SCRA,NOCC(8+I),NBFIRR(I),IUHF)
#endif
                  CALL SQUEZ2(SCRA,DENS(LDIM1+ITRIOF(I)),NBFIRR(I))
  820          CONTINUE
           ENDIF
      ENDIF
C 
#ifdef _DEBUG_LVLM
      Write(6,*) 
      Write(6, "(a)") " The checksums of DENS matrices"
      call checksum("DENSA", DENS(1), LDIM1)
      call checksum("DENSB", DENS(LDIM1), LDIM1)
      Write(6,*)
#endif
      IF (ROHF) THEN
         IF (IFLAGS(1) .GE. 10) WRITE(LUOUT,2050)
 2050    FORMAT(' @INITGES-I, Writing out MOVECTOR record. ')
C
         CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,SCR1)
         CALL   ZERO(SCRN,NBAS*NBAS)
         DO 830 I=1,NIRREP
            CALL GETBLK(SCR1,SCRN(ISQROF(I)),NBFIRR(I),NBAS,IREPS(I))
  830    CONTINUE

         cALL PUTREC(20,'JOBARC','MOVECTOR',ISQRLN(NIRREP+1)*IINTFP,
     &               SCRN)
         ROHFMO = .TRUE.
      ENDIF
C
      RETURN
      END

