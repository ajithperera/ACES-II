      SUBROUTINE GIAONE(MAXA,MAXB,ISTEPA,ISTEPB,NUCAB,IOFF,IOFF0,ITYPE,
     *                  IAB0,WORD,IPRINT,CXYZ,SIGN)
C
C  CALCULATES EXPANSION COEFFICIENTS REQUIRED FOR GIAO CALCULATIONS
C
CEND
C
C  SEP/91 AND NOV/91 JG                       
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      CHARACTER WORD*4
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
#include "lwork.par"
C
      DIMENSION CXYZ(MXAOSQ) 
C
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /TWOVEC/ DIFPAX(MXAOSQ), DIFPAY(MXAOSQ), DIFPAZ(MXAOSQ),
     *                DIFPBX(MXAOSQ), DIFPBY(MXAOSQ), DIFPBZ(MXAOSQ),
     *                TEXP1(MXAOSQ),  TEXP2(MXAOSQ),  HEXPPI(MXAOSQ)
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
      IF (ITYPE .EQ. 1) THEN
       MAX1 = MAXA
       MAX2 = MAXB
       ISTEP1 = ISTEPA
       ISTEP2 = ISTEPB
       IADDTE = 0
      ELSE
       MAX1 = MAXB
       MAX2 = MAXA
       ISTEP1 = ISTEPB
       ISTEP2 = ISTEPA
       IADDTE = MXAOSQ
      END IF
C
      IAB0P = IAB0 + 1
      IADD = - IOFF + IOFF0 + NUCAB
      IADDNP = - NUCAB - NUCAB
      NUCABD = IAB0P*NUCAB
C
      IADR00 = IOFF
      DO 200 I1 = 0, MAX1
       IADR0 = IADR00
       DO 210 I2 = 0, MAX2
        I12 = I1 + I2
        MINT = IBTAND(I12 + 1,IAB0)
        IADR = IADR0 + MINT*NUCAB
        INEXT = IADR + IADD
        IPREV = INEXT + IADDNP
        ISAME = INEXT - NUCAB
        DO 220 IT = MINT, I12 + 1, IAB0P
         IF(IT.EQ.I12) THEN
          IF(IT.NE.0) THEN
          DO 120 I=1,NUCAB
           WORK2(IADR+I)=CXYZ(I)*WORK2(ISAME+I)+HEXPPI(I)*
     &                    WORK2(IPREV+I)
120       CONTINUE
          ELSE
          DO 121 I=1,NUCAB
           WORK2(IADR+I)=CXYZ(I)*WORK2(ISAME+I)
121       CONTINUE
          ENDIF
         ELSE IF(IT.EQ.0) THEN
          DO 122 I=1,NUCAB
           WORK2(IADR+I)=CXYZ(I)*WORK2(ISAME+I)+SIGN*FLOAT(IT+1)*
     &                   WORK2(INEXT+I)
122       CONTINUE
         ELSE IF(IT.LE.(I12-1)) THEN
         DO 230 I = 1, NUCAB
          WORK2(IADR + I) = HEXPPI(I)*WORK2(IPREV+I)+CXYZ(I)*
     &                      WORK2(ISAME+I)+SIGN*FLOAT(IT+1)*
     &                      WORK2(INEXT+I)
230      CONTINUE
         ELSE IF(IT.EQ.(I12+1)) THEN
          DO 250 I=1,NUCAB
          WORK2(IADR + I) = HEXPPI(I)*WORK2(IPREV+I)
250      CONTINUE
         ELSE
          CALL ERREX
         ENDIF
         IADR = IADR + NUCABD
         INEXT = INEXT + NUCABD
         IPREV = IPREV + NUCABD
         ISAME = ISAME + NUCABD
  220   CONTINUE
        IADR0 = IADR0 + ISTEP2
  210  CONTINUE
       IADR00 = IADR00 + ISTEP1
  200 CONTINUE
      IF (IPRINT .LT. 10) RETURN
C
C  PRINT SECTION 
C
      WRITE (LUPRI, 2000)
      WRITE (LUPRI, 2010) MAX1, MAX2
      WRITE (LUPRI, 2020) ISTEPA, ISTEPB
      WRITE (LUPRI, 2030) NUCAB
      WRITE (LUPRI, 2040) IOFF0
      WRITE (LUPRI, 2050) IOFF
      WRITE (LUPRI, 2060) ITYPE
      IF (IPRINT .LT. 20) RETURN
      IADR0 = IOFF
      DO 1000 IA = 0, MAXA
       IADR = IADR0
       DO 1100 IB = 0, MAXB
        IADRT = IADR
        DO 1200 IT = 0, IA + IB + 1
         IODD = IBTAND(IA + IB + 1 - IT,IAB0)
         IF (IODD .EQ. 0) THEN
          WRITE (LUPRI, 2070) WORD, IA, IB, IT
          WRITE (LUPRI, 2080) (WORK2(IADRT + I), I = 1, NUCAB)
         END IF
         IADRT = IADRT + NUCAB
 1200   CONTINUE
        IADR = IADR + ISTEPB
 1100  CONTINUE
       IADR0 = IADR0 + ISTEPA
 1000 CONTINUE
      RETURN
 2000 FORMAT (//,'  <<<<<<<<<< SUBROUTINE GIAONE >>>>>>>>>> ',/)
 2010 FORMAT ('  MAXA/B:   ',2I7)
 2020 FORMAT ('  ISTEPA/B: ',2I7)
 2030 FORMAT ('  NUCAB:    ',I7)
 2040 FORMAT ('  IOFF0:    ',I7)
 2050 FORMAT ('  IOFF:     ',I7)
 2060 FORMAT ('  ITYPE:    ',I7)
 2070 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
 2080 FORMAT (1X,6F12.8)
      END
