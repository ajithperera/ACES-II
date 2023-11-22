
      SUBROUTINE GEN_CORE_HAM_OCC(SCR1, SCRN, ONEH, EVAL, SCRA,
     &                            SCR2, SCRB, XFORM, SCRTMP,
     &                            NOCC, LDIM1, LDIM2, IUHF,
     &                            NBAS)

       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
#include "molcas.com"
#include "symm2.com"
#include "machsp.com"
C
      DIMENSION SCR1(NBAS*NBAS),SCRN(NBAS*NBAS),ONEH(LDIM1),
     &           XFORM(LDIM1),SCR2(LDIM2),SCRTMP(LDIM2),
     &           SCRA(LDIM2),SCRB(LDIM2), NOCC(16),
     &           EVAL((IUHF+1)*NBAS)
C
      INDX2(I,J,N)=I+(J-1)*N
C
      CALL ZERO(SCR1,NBAS*NBAS)
      CALL ZERO(SCRN,NBAS*NBAS)
C
      ICNT = 0
      DO I=1,NIRREP
         IF(NBFIRR(I).NE. 0) THEN
C
           CALL EXPND2(ONEH(ITRIOF(I)),SCR2,NBFIRR(I))
           CALL EXPND2(XFORM(ITRIOF(I)),SCRA,NBFIRR(I))
           CALL TRANSP(SCRA,SCRB,NBFIRR(I),NBFIRR(I))
           CALL MXM(SCRB,NBFIRR(I),SCR2,NBFIRR(I),SCRTMP,NBFIRR(I))
           CALL MXM(SCRTMP,NBFIRR(I),SCRA,NBFIRR(I),SCR2,NBFIRR(I))
           CALL ZERO(SCR1,NBAS*NBAS)
           CALL EIG(SCR2,SCR1,NBFIRR(I),NBFIRR(I),0)
           DO  J=1,NBFIRR(I)
               ICNT=ICNT+1
               EVAL(ICNT)=SCR2(INDX2(J,J,NBFIRR(I)))
           ENDDO
         ENDIF
C
      ENDDO
C
      IF (IUHF.EQ.1) THEN
         CALL DCOPY(NBAS,EVAL(1),1,EVAL(NBAS+1),1)
      ENDIF
C
      CALL OCCUPY(NIRREP,NBFIRR,NBAS,EVAL,SCR1,IUHF)
C
      WRITE(6,5000)
      WRITE(6,5001)(NOCC(I),I=1,NIRREP)
      WRITE(6,5002)(NOCC(8+I),I=1,NIRREP)
      WRITE(6,5010)
 5000 FORMAT(T3,'Occupancies from core Hamiltonian:',/)
 5001 FORMAT(T8,'   Alpha population by irrep: ',8(I3,2X))
 5002 FORMAT(T8,'    Beta population by irrep: ',8(I3,2X))
 5010 FORMAT(/)
C
      RETURN
      END

