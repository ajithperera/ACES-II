      SUBROUTINE EXPORT_AODENS(WORK, DENS, SCR2, LDIM1, LDIM2, MAXCOR,
     &                         NBASX, NBAS, IUHF)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

#include "machsp.com"
#include "symm2.com"
C
      DIMENSION WORK(MAXCOR), DENS((IUHF+1)*LDIM1), SCR2(LDIM2)
C
      COMMON /FLAGS/  IFLAGS(100)
C
C  Read AO density from AODENS file.
C
      OPEN(UNIT=61,FILE='AODENS',STATUS='OLD',ACCESS='SEQUENTIAL',
     &     FORM='FORMATTED')
C
      I000 = 1
      I010 = I000 + NBASX*NBASX
      I020 = 1010 + NBASX*NBASX
      I030 = I020 + NBASX*NBAS 
      IEND = I030 + NBASX*NBASX
C
      CALL GETAODENS(WORK(I000), 61, NBASX)
      IF (IUHF .EQ. 1)  CALL GETAODENS(WORK(I010), 61, NBASX)
        
      CLOSE(61,STATUS='KEEP')
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The AO density matrix in EXPORT_GUESS_MO_ORBS"
      call output(work(i000), 1, nbasx, 1, nbasx, nbasx, nbasx, 1)
      if (iuhf .eq. 1) call output(work(i010), 1, nbasx, 1, nbasx,
     &                             nbasx, nbasx, 1)
#endif
C
      CALL GETREC(20,'JOBARC','ZMAT2CMP',NBASX*NBAS*IINTFP,
     &             WORK(I020))
      CALL XGEMM('N','N',NBAS,NBASX,NBASX,1.0D+00,WORK(I020),
     &             NBAS,WORK(I000),NBASX,0.0D+00,WORK(I030),NBAS)
      CALL XGEMM('N','T',NBAS,NBAS,NBASX,1.0D+00,WORK(I020),NBAS,
     &            WORK(I030),NBAS,0.0D+00,WORK(I000),NBAS)
      IF (IUHF .EQ. 1) THEN
         CALL XGEMM('N','N',NBAS,NBASX,NBASX,1.0D+00,WORK(I020),
     &               NBAS,WORK(I010),NBASX,0.0D+00,WORK(I030),NBAS)
         CALL XGEMM('N','T',NBAS,NBAS,NBASX,1.0D+00,WORK(I020),NBAS,
     &               WORK(I030),NBAS,0.0D+00,WORK(I010),NBAS)
      ENDIF
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The sym ada. density matrix in EXPORT_AODENS"
      call output(work(i000), 1, nbas, 1, nbas, nbas, nbas, 1)
      if (iuhf. eq. 1) call output(work(i010), 1, nbas, 1, nbas, nbas, 
     &                             nbas, 1)
#endif

      DO I = 1, NIRREP
         CALL GETBLK(WORK(I000), SCR2, NBFIRR(I), NBAS, IREPS(I))
         CALL SQUEZ2(SCR2, DENS(ITRIOF(I)), NBFIRR(I))
C
#ifdef _DEBUG_LVL0
        Write(6,*) "The symmetry packed density in EXPORT_AODENS"
        call output(scr2, 1, nbfirr(i), 1, nbfirr(i), nbfirr(i),
     &              nbfirr(i), 1)
#endif
C
         IF (IUHF .EQ. 1) THEN
            CALL GETBLK(WORK(I010), SCR2, NBFIRR(I), NBAS, IREPS(I))
#ifdef _DEBUG_LVL0
        Write(6,*) "The symmetry packed density in EXPORT_AODENS"
        call output(scr2, 1, nbfirr(i), 1, nbfirr(i), nbfirr(i),
     &              nbfirr(i), 1)
#endif
            CALL SQUEZ2(SCR2, DENS(ITRIOF(I)+LDIM1), NBFIRR(I))
         ENDIF
      ENDDO
c
      RETURN
      END

