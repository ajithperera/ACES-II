      SUBROUTINE CNTOUT(TWOIN,TWOOUT,WORK1,WORK2,LWORK1,LWORK2,
     &                  ISTRCF,IOFFVC,NSET1,NSET2,
     &                  NPRIM1,NPRIM2,NCONT1,NCONT2,
     &                  NINT34,TPRI12,TCON12)
C
C
C     Purpose: Transformation of two outermost indices
C              Index 1 is outermost index
C              Index 2 is next outermost index
C
C     Note: After transformation the order of two outermost indices
C           is reversed
C
C     In:  TWOIN(NPRIM1*NPRIM2*NINT34)
C
C     Out: TWOOUT(NCONT2*NCONT1*NINT34)
C
C     Scratch: WORK1(NCONT1*NPRIM2*NINT34)
C
C              WORK2(NPRIM1*NPRIM2*NINT34)              (TPRI12)
C              WORK2(NCONT1*NPRIM2*NINT34)         (.NOT.TPRI12)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "baslims.par"
      COMMON /CCFCOM/ CONT1 (MXCONT*MXAOVC), CONT2 (MXCONT*MXAOVC),
     &                CONT3 (MXCONT*MXAOVC), CONT4 (MXCONT*MXAOVC),
     &                CONTT1(MXCONT*MXAOVC), CONTT2(MXCONT*MXAOVC),
     &                CONTT3(MXCONT*MXAOVC), CONTT4(MXCONT*MXAOVC),
     &                NUC1X (MXAOVC),        NUC2X (MXAOVC),
     &                NUC3X (MXAOVC),        NUC4X (MXAOVC),
     &                NRC1X (MXAOVC),        NRC2X (MXAOVC),
     &                NRC3X (MXAOVC),        NRC4X (MXAOVC)
C
      LOGICAL TPRI12, TCON12
      DIMENSION TWOIN(NPRIM1*NPRIM2*NINT34),
     &          TWOOUT(NCONT1*NCONT2*NINT34),
     &          WORK1(LWORK1),WORK2(LWORK2)
C
      DATA AZERO,ONE /0.D0,1.D0/
C
      NPR234 = NPRIM2*NINT34
      NCT134 = NCONT1*NINT34
C
C TRANSFORM FIRST INDEX
C
      IF (TPRI12) THEN
       DO 100 I = 1, NPRIM1
        IJ = (I - 1)*NPRIM1*NINT34
        IADR = ((I-1)*I)/2*NINT34
        DO 110 K = 1, I*NINT34
          WORK2(IJ + K) = TWOIN(IADR + K)
  110    CONTINUE
  100  CONTINUE
       IADR = NINT34
       DO 200 I = 2, NPRIM1
        DO 210 J = 1, I - 1
         JI = ((J - 1)*NPRIM1 + I - 1)*NINT34
         DO 220 K = 1, NINT34
          WORK2(JI + K) = TWOIN(IADR + K)
  220    CONTINUE
         IADR = IADR + NINT34
  210   CONTINUE
        IADR = IADR + NINT34
  200  CONTINUE
       IWORK1 = 1
       IWORK2 = 1
       ICONT  = ISTRCF
       DO 300 I = 1, NSET1
        NUC1 = NUC2X(IOFFVC + I)
        NRC1 = NRC2X(IOFFVC + I)
        CALL XGEMM('N','N',NPR234,NRC1,NUC1,ONE,WORK2(IWORK2),NPR234,
     &             CONTT2(ICONT),NUC1,AZERO,WORK1(IWORK1),NPR234)
        IWORK1 = IWORK1 + NRC1*NPR234
        IWORK2 = IWORK2 + NUC1*NPR234
        ICONT  = ICONT  + NRC1*NUC1
  300  CONTINUE
      ELSE
       IWORK1 = 1
       IWORK2 = 1
       ICONT  = ISTRCF
       DO 400 I = 1, NSET1
        NUC1 = NUC2X(IOFFVC + I)
        NRC1 = NRC2X(IOFFVC + I)
        CALL XGEMM('N','N',NPR234,NRC1,NUC1,ONE,TWOIN(IWORK2),NPR234,
     &             CONTT2(ICONT),NUC1,AZERO,WORK1(IWORK1),NPR234)
        IWORK1 = IWORK1 + NRC1*NPR234
        IWORK2 = IWORK2 + NUC1*NPR234
        ICONT  = ICONT  + NRC1*NUC1
  400  CONTINUE
      END IF
C
C CHANGE ORDER OF FIRST AND SECOND INDICES
C
      IADR1  = 0
      IADR20 = 0
      DO 500 I = 1, NCONT1
       IADR2 = IADR20
       DO 510 J = 1, NPRIM2
        DO 520 K = 1, NINT34
         WORK2(IADR2+K) = WORK1(IADR1 + K)
  520   CONTINUE
        IADR1 = IADR1 + NINT34
        IADR2 = IADR2 + NCT134
  510  CONTINUE
       IADR20 = IADR20 + NINT34
  500 CONTINUE
C
C TRANSFORM SECOND INDEX
C
      IF (TCON12) THEN
       IWORK1 = 1
       IWORK2 = 1
       ICONT  = ISTRCF
       DO 600 I = 1, NSET2
        NUC2 = NUC1X(IOFFVC + I)
        NRC2 = NRC1X(IOFFVC + I)
        CALL XGEMM('N','N',NCT134,NRC2,NUC2,ONE,WORK2(IWORK2),NCT134,
     &             CONTT1(ICONT),NUC2,AZERO,WORK1(IWORK1),NCT134)
        IWORK1 = IWORK1 + NRC2*NCT134
        IWORK2 = IWORK2 + NUC2*NCT134
        ICONT  = ICONT  + NRC2*NUC2
  600  CONTINUE
       DO 700 I = 1, NCONT1
        IADR1 = ((I - 1)*I)/2*NINT34
        IADR2 = (I - 1)*NCT134
        DO 710 J = 1, I*NINT34
         TWOOUT(IADR1 + J) = WORK1(IADR2 + J)
  710   CONTINUE
  700  CONTINUE
      ELSE
       IWORK1 = 1
       IWORK2 = 1
       ICONT  = ISTRCF
       DO 800 I = 1, NSET2
        NUC2 = NUC1X(IOFFVC + I)
        NRC2 = NRC1X(IOFFVC + I)
        CALL XGEMM('N','N',NCT134,NRC2,NUC2,ONE,WORK2(IWORK2),NCT134,
     &             CONTT1(ICONT),NUC2,AZERO,TWOOUT(IWORK1),NCT134)
        IWORK1 = IWORK1 + NRC2*NCT134
        IWORK2 = IWORK2 + NUC2*NCT134
        ICONT  = ICONT  + NRC2*NUC2
  800  CONTINUE
      END IF
      RETURN
      END
