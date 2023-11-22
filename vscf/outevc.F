
      SUBROUTINE OUTEVC(EVAL,EVEC,NBAS,LDIM1,IUHF,AOFIL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include "maxbasfn.par"
      PARAMETER(LUOUT=6)
      PARAMETER(LUINT=10)
C
      LOGICAL AOFIL
      INTEGER NBAS
      CHARACTER*4 ATMLBL(MAXBASFN),ANGLBL(MAXBASFN),junk
      CHARACTER*5 SPIN(2)
      CHARACTER*80 FNAME
      DIMENSION EVAL((IUHF+1)*NBAS),EVEC((IUHF+1)*LDIM1)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
#include "symm2.com"
C
      DATA SPIN /'ALPHA',' BETA'/
C
      INDX2(I,J,N)=I+(J-1)*N
C
C  PICK UP THE LABELS FROM THE INTEGRAL FILE.
C
      IF(AOFIL) THEN
       CALL GFNAME('IIII    ',FNAME,ILENGTH)
       OPEN(LUINT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',ACCESS=
     &     'SEQUENTIAL')
      ENDIF
      CALL LOCATE(LUINT,'LABBASIS')
       DO I=1,NBAS
          READ(LUINT) J,ATMLBL(J),ANGLBL(J),IJUNK
       END DO
C
      ICNT0=0
      DO 100 ISPIN=1,(IUHF+1)
        ICNT1=0
        WRITE(LUOUT,5001)
 5001   FORMAT(/,T3,'ORBITAL EIGENVECTORS',/)
        DO 105 I=1,NIRREP
          WRITE(LUOUT,5010)I,SPIN(ISPIN)
 5010     FORMAT(T3,'SYMMETRY BLOCK ',I1,' (',A5,')',/)
          NBY4=NBFIRR(I)/4
          DO 110 J=1,NBY4
            ICNT0=4*J+ICNT1
            WRITE(LUOUT,5019)(ICNT0-4+K,K=1,4)
 5019       FORMAT(T25,2X,'MO #',I3,3(7X,'MO #',I3),/)
            WRITE(LUOUT,5020)(EVAL((ISPIN-1)*NBAS+IRPOFF(I)+
     &                             (J-1)*4+K),K=1,4)
 5020       FORMAT(T5,'BASIS/ORB E',9X,F10.5,4X,F10.5,4X,F10.5,
     &             4X,F10.5)
            WRITE(LUOUT,5025)
 5025       FORMAT(T25,'----------',4X,'----------',4X,
     &                 '----------',4X,'----------')
            DO 115 K=1,NBFIRR(I)
              WRITE(LUOUT,5030)K+IRPOFF(I),ATMLBL(IRPOFF(I)+K),
     &                         ANGLBL(IRPOFF(I)+K),
     &                         (EVEC((ISPIN-1)*LDIM1+ISQROF(I)-1+
     &                               INDX2(K,(J-1)*4+L,NBFIRR(I))),
     &                          L=1,4)
 5030         FORMAT(T3,I3,';',1X,A4,2X,A4,7X,F10.5,4X,F10.5,4X,
     &               F10.5,4X,F10.5)
  115       CONTINUE
            WRITE(LUOUT,5035)
 5035       FORMAT(T3,/)
  110     CONTINUE
C
C  Now take care of the remaining eigenvectors, if we had a total
C  in this irrep which was not divisible by 4.
C
          IF(NBY4.EQ.0) ICNT0=ICNT1
          IF(4*NBY4.NE.NBFIRR(I)) THEN
            ICNT=NBFIRR(I)-4*NBY4
            IF(ICNT.EQ.1) THEN
              WRITE(LUOUT,5050)ICNT0+1
              WRITE(LUOUT,5021)(EVAL((ISPIN-1)*NBAS+IRPOFF(I)+
     &                               NBY4*4+K),K=1,1)
              WRITE(LUOUT,5026)
              DO 116 K=1,NBFIRR(I)
                WRITE(LUOUT,5031)K+IRPOFF(I),ATMLBL(IRPOFF(I)+K),
     &                           ANGLBL(IRPOFF(I)+K),
     &                           (EVEC((ISPIN-1)*LDIM1+ISQROF(I)-1+
     &                                 INDX2(K,NBY4*4+L,NBFIRR(I))),
     &                            L=1,1)
  116         CONTINUE
              ICNT1=ICNT0+1
            ELSEIF(ICNT.EQ.2) THEN
              WRITE(LUOUT,5051)(ICNT0+K,K=1,2)
              WRITE(LUOUT,5022)(EVAL((ISPIN-1)*NBAS+IRPOFF(I)+
     &                               NBY4*4+K),K=1,2)
              WRITE(LUOUT,5027)
              DO 117 K=1,NBFIRR(I)
                WRITE(LUOUT,5032)K+IRPOFF(I),ATMLBL(IRPOFF(I)+K),
     &                           ANGLBL(IRPOFF(I)+K),
     &                           (EVEC((ISPIN-1)*LDIM1+ISQROF(I)-1+
     &                                 INDX2(K,NBY4*4+L,NBFIRR(I))),
     &                            L=1,2)
  117         CONTINUE
              ICNT1=ICNT0+2
            ELSE
              WRITE(LUOUT,5052)(ICNT0+K,K=1,3)
              WRITE(LUOUT,5023)(EVAL((ISPIN-1)*NBAS+IRPOFF(I)+
     &                               NBY4*4+K),K=1,3)
              WRITE(LUOUT,5028)
              DO 118 K=1,NBFIRR(I)
                WRITE(LUOUT,5033)K+IRPOFF(I),ATMLBL(IRPOFF(I)+K),
     &                           ANGLBL(IRPOFF(I)+K),
     &                           (EVEC((ISPIN-1)*LDIM1+ISQROF(I)-1+
     &                                 INDX2(K,NBY4*4+L,NBFIRR(I))),
     &                            L=1,3)
  118         CONTINUE
              ICNT1=ICNT0+3
            ENDIF
 5021       FORMAT(T5,'BASIS\\ORB E',9X,F10.5,1(4X,F10.5))
 5022       FORMAT(T5,'BASIS\\ORB E',9X,F10.5,2(4X,F10.5))
 5023       FORMAT(T5,'BASIS\\ORB E',9X,F10.5,3(4X,F10.5))
 5026       FORMAT(T25,'----------')
 5027       FORMAT(T25,'----------',4X,'----------')
 5028       FORMAT(T25,'----------',4X,'----------',4X,'----------')
 5031       FORMAT(T3,I3,';',1X,A4,2X,A4,7X,F10.5,1(4X,F10.5))
 5032       FORMAT(T3,I3,';',1X,A4,2X,A4,7X,F10.5,2(4X,F10.5))
 5033       FORMAT(T3,I3,';',1X,A4,2X,A4,7X,F10.5,3(4X,F10.5))
 5050       FORMAT(T25,2X,'MO #',I3,/)
 5051       FORMAT(T25,2X,'MO #',I3,7X,'MO #',I3,/)
 5052       FORMAT(T25,2X,'MO #',I3,7X,'MO #',I3,7X,'MO #',I3,/)
            WRITE(LUOUT,5035)
          ELSE
            ICNT1=ICNT0
          ENDIF
  105   CONTINUE
  100 CONTINUE
C
      RETURN
      END
