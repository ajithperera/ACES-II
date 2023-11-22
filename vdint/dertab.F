      SUBROUTINE DERTAB(MAXA,MAXB,ISTEPA,ISTEPB,NUCAB,IOFF,IOFF0,
     *                  IAB0,WORD,IPRINT)
C
C     TUH 84
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      CHARACTER WORD*4
#include "lwork.par"
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /TWOVEC/ DIFPAX(MXAOSQ), DIFPAY(MXAOSQ), DIFPAZ(MXAOSQ),
     *                DIFPBX(MXAOSQ), DIFPBY(MXAOSQ), DIFPBZ(MXAOSQ),
     *                TEXP1(MXAOSQ),  TEXP2(MXAOSQ),  HEXPPI(MXAOSQ)
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
      IAB0P = IAB0 + 1
      NUCADD = IAB0P*NUCAB
      IADD0 = - IOFF + IOFF0
      IADDMM = IADD0 - ISTEPA - ISTEPB
      IADDMP = IADD0 - ISTEPA + ISTEPB
      IADDPM = IADD0 + ISTEPA - ISTEPB
      IADDPP = IADD0 + ISTEPA + ISTEPB
C
C     ************************************
C     ********** IA = 0, IB = 0 **********
C     ************************************
C
      ISTRAB = IOFF
      IF (IAB0 .EQ. 0) THEN
         IADR0 = ISTRAB
         IADR1 = IADR0 + NUCADD
         IADR2 = IADR1 + NUCADD
         IODPP0 = IADR0 + IADDPP
         IODPP1 = IADR1 + IADDPP
         IODPP2 = IADR2 + IADDPP
         DO 100 I = 1, NUCAB
            FACPP = TEXP1(I)*TEXP2(I)
            WORK2(IADR0 + I) = FACPP*WORK2(IODPP0 + I)
            WORK2(IADR1 + I) = FACPP*WORK2(IODPP1 + I)
            WORK2(IADR2 + I) = FACPP*WORK2(IODPP2 + I)
  100    CONTINUE
      ELSE
         IADR0 = ISTRAB
         IADR1 = IADR0 + NUCADD
         IODPP0 = IADR0 + IADDPP
         IODPP1 = IADR1 + IADDPP
         DO 105 I = 1, NUCAB
            FACPP = TEXP1(I)*TEXP2(I)
            WORK2(IADR0 + I) = FACPP*WORK2(IODPP0 + I)
            WORK2(IADR1 + I) = FACPP*WORK2(IODPP1 + I)
  105    CONTINUE
      END IF
C
C     ************************************
C     ********** IA = 0, IB > 0 **********
C     ************************************
C
      ISTRAB = ISTRAB + ISTEPB
      DO 200 IB = 1, MAXB
         FACPM = - FLOAT(IB)
C
C        ***** IT = 0, IB *****
C
         MINT = IBTAND(IB,IAB0)
         IADR0 = ISTRAB + MINT*NUCAB
         DO 210 IT = MINT, IB, IAB0P
            IODPP0 = IADR0 + IADDPP
            IODPM0 = IADR0 + IADDPM
            DO 220 I = 1, NUCAB
               TEXPA = TEXP1(I)
               TEXPB = TEXP2(I)
               WORK2(IADR0 + I) = TEXPA*(TEXPB*WORK2(IODPP0 + I)
     *                                + FACPM*WORK2(IODPM0 + I))
  220       CONTINUE
            IADR0 = IADR0 + NUCADD
  210    CONTINUE
C
C        ***** IT = IB + 1, IB + 2 *****
C
         IF (IAB0 .EQ. 0) THEN
            IADR1 = IADR0 + NUCADD
            IODPP0 = IADR0 + IADDPP
            IODPP1 = IADR1 + IADDPP
            DO 230 I = 1, NUCAB
               FACPP = TEXP1(I)*TEXP2(I)
               WORK2(IADR0 + I) = FACPP*WORK2(IODPP0 + I)
               WORK2(IADR1 + I) = FACPP*WORK2(IODPP1 + I)
  230       CONTINUE
         ELSE
            IODPP0 = IADR0 + IADDPP
            DO 235 I = 1, NUCAB
               WORK2(IADR0 + I) = TEXP1(I)*TEXP2(I)*WORK2(IODPP0 + I)
  235       CONTINUE
         END IF
         ISTRAB = ISTRAB + ISTEPB
  200 CONTINUE
C
C     ****************************
C     ********** IA > 0 **********
C     ****************************
C
      ISTRTA = IOFF + ISTEPA
      DO 300 IA = 1, MAXA
         FACMP = - FLOAT(IA)
         ISTRAB = ISTRTA
C
C        ****************************
C        ********** IB = 0 **********
C        ****************************
C
C        ***** IT = 0, IA *****
C
         MINT = IBTAND(IA,IAB0)
         IADR0 = ISTRAB + MINT*NUCAB
         DO 310 IT = MINT, IA, IAB0P
            IODPP0 = IADR0 + IADDPP
            IODMP0 = IADR0 + IADDMP
            DO 320 I = 1, NUCAB
               TEXPA = TEXP1(I)
               TEXPB = TEXP2(I)
               WORK2(IADR0 + I) = TEXPB*(TEXPA*WORK2(IODPP0 + I)
     *                                + FACMP*WORK2(IODMP0 + I))
  320       CONTINUE
            IADR0 = IADR0 + NUCADD
  310    CONTINUE
C
C        ***** IT = IA + 1, IA + 2 *****
C
         IF (IAB0 .EQ. 0) THEN
            IADR1 = IADR0 + NUCADD
            IODPP0 = IADR0 + IADDPP
            IODPP1 = IADR1 + IADDPP
            DO 330 I = 1, NUCAB
               FACPP = TEXP1(I)*TEXP2(I)
               WORK2(IADR0 + I) = FACPP*WORK2(IODPP0 + I)
               WORK2(IADR1 + I) = FACPP*WORK2(IODPP1 + I)
  330       CONTINUE
         ELSE
            IODPP0 = IADR0 + IADDPP
            DO 335 I = 1, NUCAB
               WORK2(IADR0 + I) = TEXP1(I)*TEXP2(I)*WORK2(IODPP0 + I)
  335       CONTINUE
         END IF
C
C        ****************************
C        ********** IB > 0 **********
C        ****************************
C
         ISTRAB = ISTRAB + ISTEPB
         DO 400 IB = 1, MAXB
            FACPM = - FLOAT(IB)
            FACMM = FLOAT(IA*IB)
            IAB = IA + IB
C
C           ***** IT =  0, IAB - 2 *****
C
            MINT = IBTAND(IAB,IAB0)
            IADR0 = ISTRAB + MINT*NUCAB
            DO 410 IT = MINT, IAB - 2, IAB0P
               IODPP0 = IADR0 + IADDPP
               IODPM0 = IADR0 + IADDPM
               IODMP0 = IADR0 + IADDMP
               IODMM0 = IADR0 + IADDMM
               DO 420 I = 1, NUCAB
                  TEXPA = TEXP1(I)
                  TEXPB = TEXP2(I)
                  WORK2(IADR0 + I) = TEXPA*TEXPB*WORK2(IODPP0 + I)
     *                             + TEXPA*FACPM*WORK2(IODPM0 + I)
     *                             + TEXPB*FACMP*WORK2(IODMP0 + I)
     *                             + FACMM*WORK2(IODMM0 + I)
  420          CONTINUE
               IADR0 = IADR0 + NUCADD
  410       CONTINUE
C
C           ***** IT = IAB - 1, IAB + 2 *****
C
            IF (IAB0 .EQ. 0) THEN
               IADR1 = IADR0 + NUCADD
               IADR2 = IADR1 + NUCADD
               IADR3 = IADR2 + NUCADD
               IODPM0 = IADR0 + IADDPM
               IODPM1 = IADR1 + IADDPM
               IODMP0 = IADR0 + IADDMP
               IODMP1 = IADR1 + IADDMP
               IODPP0 = IADR0 + IADDPP
               IODPP1 = IADR1 + IADDPP
               IODPP2 = IADR2 + IADDPP
               IODPP3 = IADR3 + IADDPP
               DO 430 I = 1, NUCAB
                  TEXPA = TEXP1(I)
                  TEXPB = TEXP2(I)
                  FACPP = TEXPA*TEXPB
                  TFACPM = TEXPA*FACPM
                  TFACMP = TEXPB*FACMP
                  WORK2(IADR0 + I) = FACPP*WORK2(IODPP0 + I)
     *                            + TFACPM*WORK2(IODPM0 + I)
     *                            + TFACMP*WORK2(IODMP0 + I)
                  WORK2(IADR1 + I) = FACPP*WORK2(IODPP1 + I)
     *                            + TFACPM*WORK2(IODPM1 + I)
     *                            + TFACMP*WORK2(IODMP1 + I)
                  WORK2(IADR2 + I) = FACPP*WORK2(IODPP2 + I)
                  WORK2(IADR3 + I) = FACPP*WORK2(IODPP3 + I)
  430          CONTINUE
            ELSE
               IADR1 = IADR0 + NUCADD
               IODPP0 = IADR0 + IADDPP
               IODPP1 = IADR1 + IADDPP
               IODPM0 = IADR0 + IADDPM
               IODMP0 = IADR0 + IADDMP
               DO 435 I = 1, NUCAB
                  TEXPA = TEXP1(I)
                  TEXPB = TEXP2(I)
                  FACPP = TEXPA*TEXPB
                  WORK2(IADR0 + I) = FACPP*WORK2(IODPP0 + I)
     *                             + TEXPA*FACPM*WORK2(IODPM0 + I)
     *                             + TEXPB*FACMP*WORK2(IODMP0 + I)
                  WORK2(IADR1 + I) = FACPP*WORK2(IODPP1 + I)
  435          CONTINUE
            END IF
            ISTRAB = ISTRAB + ISTEPB
  400    CONTINUE
         ISTRTA = ISTRTA + ISTEPA
  300 CONTINUE
      IF (IPRINT .LT. 10) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      WRITE (LUPRI, 2000)
      WRITE (LUPRI, 2010) MAXA, MAXB
      WRITE (LUPRI, 2020) ISTEPA, ISTEPB
      WRITE (LUPRI, 2030) NUCAB
      WRITE (LUPRI, 2040) IOFF0
      WRITE (LUPRI, 2050) IOFF
      IF (IPRINT .LT. 20) RETURN
      ISTRAB = IOFF
      DO 1000 IA = 0, MAXA
         IADR = ISTRAB
         DO 1100 IB = 0, MAXB
            IADRT = IADR
            DO 1200 IT = 0, IA + IB + 2
               IODD = IBTAND(IA + IB + 2 - IT,IAB0)
               IF (IODD .EQ. 0) THEN
                  WRITE (LUPRI, 2070) WORD, IA, IB, IT
                  WRITE (LUPRI, 2080) (WORK2(IADRT + I), I = 1, NUCAB)
               END IF
               IADRT = IADRT + NUCAB
 1200       CONTINUE
            IADR = IADR + ISTEPB
 1100    CONTINUE
         ISTRAB = ISTRAB + ISTEPA
 1000 CONTINUE
      RETURN
 2000 FORMAT (//,'  <<<<<<<<<< SUBROUTINE DERTAB >>>>>>>>>> ',/)
 2010 FORMAT ('  MAXA/B:   ',2I7)
 2020 FORMAT ('  ISTEPA/B: ',2I7)
 2030 FORMAT ('  NUCAB:    ',I7)
 2040 FORMAT ('  IOFF0:    ',I7)
 2050 FORMAT ('  IOFF:     ',I7)
 2070 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
 2080 FORMAT (1X,6F12.8)
      END
