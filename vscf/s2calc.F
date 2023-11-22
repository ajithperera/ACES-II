      SUBROUTINE S2CALC(EVEC,SOVRLP,SCR1,SCR2,LDIM1,LDIM2,NBAS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION EVEC(2*LDIM1),SOVRLP(NBAS*NBAS)
      DIMENSION SCR1(LDIM2),SCR2(LDIM2)
      DIMENSION NOCC(16)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /POPUL/ NOCC
#include "symm2.com"
C
      DATA ONE /1.0/
      DATA TWO /2.0/
      DATA ZILCH /0.0/
C
      INDX2(I,J,N)=I+(J-1)*N
C
C  Form overlap of alpha and beta density matrices.
C
      NALPHA=0
      NBETA=0
      DO 50 I=1,NIRREP
        NALPHA=NALPHA+NOCC(I)
        NBETA=NBETA+NOCC(8+I)
   50 CONTINUE
      SAVG=(NALPHA-NBETA)/TWO
      S2=SAVG*(SAVG+ONE)
      CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS*NBAS*IINTFP,SOVRLP)
      DO I=1,NIRREP
         IF ((NOCC(I).NE.0).AND.(NOCC(8+I).NE.0)) THEN
            CALL GETBLK(SOVRLP,SCR1,NBFIRR(I),NBAS,IREPS(I))
            CALL XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),
     &                 ONE,  SCR1,           NBFIRR(I),
     &                       EVEC(ISQROF(I)),NBFIRR(I),
     &                 ZILCH,SCR2,           NBFIRR(I))
            CALL XGEMM('T','N',NOCC(8+I),NOCC(I),NBFIRR(I),
     &                 ONE,  EVEC(LDIM1+ISQROF(I)),NBFIRR(I),
     &                       SCR2,                 NBFIRR(I),
     &                 ZILCH,SCR1,                 NOCC(8+I))
            DO J=1,NOCC(I)
               S2=S2-SDOT(NOCC(8+I),SCR1(INDX2(1,J,NOCC(8+I))),1,
     &                              SCR1(INDX2(1,J,NOCC(8+I))),1)
            END DO
         END IF
      END DO
      S2=S2+NBETA
C
      AMULT=SQRT(1+4*S2)
      WRITE(LUOUT,5000)AMULT,S2
 5000 FORMAT(/,T3,'     The average multiplicity is ',F12.7,/,
     &         T3,'The expectation value of S**2 is ',F12.7,/)
      CALL PUTREC(20,'JOBARC','S2SCF   ',IINTFP,S2)
      RETURN
      END
