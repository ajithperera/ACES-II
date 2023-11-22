      SUBROUTINE OSTAU35(T2,NSIZTAR,IRREP,T1A,NSIZTA,IRREPA,POPA,
     $     T1B,NSIZTB,IRREPB,POPB,FACT,FACT0)
C
C     THIS ROUTINE CALCULATE THE OUTER PRODUCT OF VECTORS T1A AND T1B
C     AS REQURIED IN OS CALCULATIONS FOR QUANTITIES OF LIST 35+X. 
C     PLEASE NOTE THAT THE TWO FIRST OR THE TWO SECOND INDECES OF T2
C     MUST BE 1 ('AA' TYPE)
C
CEND
C
CPROGRAMED BY PS SEPT/93
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POPA,POPB,DIRPRD
      DIMENSION T2(NSIZTAR),T1A(NSIZTA),T1B(NSIZTB),POPA(8),POPB(8)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPC(255),IRREPD(255),
     &                DIRPRD(8,8)
      IOFF=1
      DO 10 IRREP2=1,IRREPB-1
         IRREP1=DIRPRD(IRREP,IRREP2)
         IOFF=IOFF+POPA(IRREP1)*POPB(IRREP2)
 10   CONTINUE
C
      DO 1 J=1,NSIZTB
         DO 2 I=1,NSIZTA
            T2(IOFF)=FACT0*T2(IOFF)+FACT*T1A(I)*T1B(J)
            IOFF=IOFF+1
 2       CONTINUE
 1    CONTINUE
C
      DO 20 IRREP2=IRREPB+1,NIRREP
         IRREP1=DIRPRD(IRREP,IRREP2)
         IOFF=IOFF+POPA(IRREP1)*POPB(IRREP2)
 20   CONTINUE
C
      IF(IOFF-1.NE.NSIZTAR) THEN
         WRITE(6,*)'@OSTAU35-E FATAL ERROR: IOFF-1.NE.NSIZTAR',
     $      IOFF-1,NSIZTAR
         CALL ERREX
      ENDIF
C
      RETURN
      END
