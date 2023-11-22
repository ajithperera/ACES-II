      SUBROUTINE HERI(JMAX,ISTEPT,ISTEPU,ISTEPV,NRTUV,
     *                IPQ0X,IPQ0Y,IPQ0Z,IPRINT,WORK1,LWORK1)
C
C    tuh fall 1984
C
C    Modified Jul 28 88 to avoid multiplying zero with
C    undetermined numbers - tuh
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      LOGICAL PQXGT0, PQYGT0, PQZGT0
      COMMON /PRIVEC/ NPRM12(MXAOSQ), NPRM34(MXAOSQ),
     *                LSTP12(MXAOSQ), LSTP34(MXAOSQ),
     *                LSTO12(MXAOSQ), LSTO34(MXAOSQ),
     *                LASTP(MXAOSQ),  LASTO(MXAOSQ)
      DIMENSION  WORK1(LWORK1)
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
#include "ninahh.par"
      COMMON /CINAHH/ INAHH(NINAHH)
#include "lwork.par"
      PARAMETER (MAXVEC = LWORK3 / 6)
      COMMON /CWORK3/ WK3LOW, DPQX(MAXVEC), DPQY(MAXVEC), DPQZ(MAXVEC),
     *                RJ000(3*MAXVEC)
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTXOR(I,J) = XOR(I,J)
C
C     R(IT,IU,IV)
C
C     R(0,0,0)  R(1,0,0)  R(2,0,0)  R(3,0,0)
C     R(0,1,0)  R(1,1,0)  R(2,1,0)    0.0
C     R(0,2,0)  R(1,2,0)    0.0       0.0
C     R(0,3,0)    0.0       0.0       0.0
C
C     R(0,0,1)  R(1,0,1)  R(2,0,1)    0.0
C     R(0,1,1)  R(1,1,1)    0.0       0.0
C     R(0,2,1)    0.0       0.0       0.0
C       0.0       0.0       0.0       0.0
C
C     R(0,0,2)  R(1,0,2)    0.0       0.0
C     R(0,1,2)    0.0       0.0       0.0
C       0.0       0.0       0.0       0.0
C       0.0       0.0       0.0       0.0
C
C     R(0,0,3)    0.0       0.0       0.0
C       0.0       0.0       0.0       0.0
C       0.0       0.0       0.0       0.0
C       0.0       0.0       0.0       0.0
C
C
      IF (JMAX .EQ. 0) THEN
         DO 50 I = 1, ISTEPT
            WORK1(IWK1HH + I) = RJ000(I)
   50    CONTINUE
         INAHH(1) = IWK1HH
         RETURN
      END IF
      PQXGT0 = IPQ0X .EQ. 0
      PQYGT0 = IPQ0Y .EQ. 0
      PQZGT0 = IPQ0Z .EQ. 0
      IF (IBTAND(JMAX,1) .EQ. 0) THEN
         ISTRTJ = IWKHHS + NRTUV
         ISTEPJ = - NRTUV
      ELSE
         ISTRTJ = IWKHHS
         ISTEPJ = NRTUV
      END IF
C
C     ***** JVAL = 1 *****
C
C     R(0,0,0)
C
      IRJMAX = ISTEPT*JMAX
      IRJ = IRJMAX - ISTEPT
      DO 100 I = 1, ISTEPT
         WORK1(ISTRTJ + I) = RJ000(IRJ + I)
  100 CONTINUE
C
C     R(1,0,0)
C
      IF (PQXGT0) THEN
         ISTART = ISTRTJ + ISTEPT
         DO 110 I = 1, ISTEPT
            WORK1(ISTART + I) = DPQX(I)*RJ000(IRJMAX + I)
  110    CONTINUE
      END IF
C
C     R(0,1,0)
C
      IF (PQYGT0) THEN
         ISTART = ISTRTJ + ISTEPU
         DO 120 I = 1, ISTEPT
            WORK1(ISTART + I) = DPQY(I)*RJ000(IRJMAX + I)
  120    CONTINUE
      END IF
C
C     R(0,0,1)
C
      IF (PQZGT0) THEN
         ISTART = ISTRTJ + ISTEPV
         DO 130 I = 1, ISTEPT
            WORK1(ISTART + I) = DPQZ(I)*RJ000(IRJMAX + I)
  130    CONTINUE
      END IF
C
C     ***** JVAL > 1 *****
C
      IF (JMAX .GT. 1) THEN
         IPQ0X1 = IPQ0X + 1
         IPQ0Y1 = IPQ0Y + 1
         IPQ0Z1 = IPQ0Z + 1
         IADDT = IPQ0X1*ISTEPT
         IADDU = IPQ0Y1*ISTEPU
         IADDV = IPQ0Z1*ISTEPV
         DO 200 JVAL = 2,JMAX
            ISTRTJ = ISTRTJ + ISTEPJ
            ISTEPJ = - ISTEPJ
C
C           R(0,0,0)
C
            IRJ = IRJ - ISTEPT
            DO 300 I = 1, ISTEPT
               WORK1(ISTRTJ + I) = RJ000(IRJ + I)
  300       CONTINUE
C
C           R(T,0,0)
C
            IF (PQXGT0) THEN
               ISTART = ISTRTJ + IADDT
               IPREV1 = ISTRTJ + ISTEPJ
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
               DO 400 I = 1, ISTEPT
                  WORK1(ISTART + I) = DPQX(I)*WORK1(IPREV1 + I)
  400          CONTINUE
               IF (JVAL .GT. 1) THEN
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
                  DO 410 IT = 2,JVAL
                     TMIN1  = FLOAT(IT - 1)
                     ISTART = ISTRTJ + IT*IADDT
                     IPREV1 = ISTART + ISTEPJ - IADDT
                     IPREV2 = IPREV1 - IADDT
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                     DO 420 I = 1, ISTEPT
                        WORK1(ISTART + I) = DPQX(I)*WORK1(IPREV1 + I)
     *                                      + TMIN1*WORK1(IPREV2 + I)
  420                CONTINUE
  410             CONTINUE
               END IF
            ELSE
               ISTART = ISTRTJ + IADDT
               IPREV2 = ISTRTJ + ISTEPJ
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
               DO 430 IT = 2, JVAL, 2
                  TMIN1 = FLOAT(IT - 1)
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                  DO 440 I = 1, ISTEPT
                     WORK1(ISTART + I) = TMIN1*WORK1(IPREV2 + I)
  440             CONTINUE
                  ISTART = ISTART + IADDT
                  IPREV2 = IPREV2 + IADDT
  430          CONTINUE
            END IF
C
C              R(T,U,0)
C
            ISTRTU = ISTRTJ + IADDU
            IADDPR = ISTEPJ - IADDU
            IF (PQYGT0) THEN
               ISTART = ISTRTU
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
               DO 500 IT = 0, JVAL - 1, IPQ0X1
                  IPREV1 = ISTART + IADDPR
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                  DO 510 I = 1, ISTEPT
                     WORK1(ISTART + I) = DPQY(I)*WORK1(IPREV1 + I)
  510             CONTINUE
                  ISTART = ISTART + IADDT
  500          CONTINUE
               IF (JVAL .GT. 1) THEN
                  ISTRTU = ISTRTU + IADDU
                  DO 520 IU = 2, JVAL
                     UMIN1 = FLOAT(IU - 1)
                     ISTART = ISTRTU
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
                     DO 530 IT = 0, JVAL - IU, IPQ0X1
                        IPREV1 = ISTART + IADDPR
                        IPREV2 = IPREV1 - IADDU
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                        DO 540 I = 1, ISTEPT
                           WORK1(ISTART + I) = DPQY(I)*WORK1(IPREV1 + I)
     *                                         + UMIN1*WORK1(IPREV2 + I)
  540                   CONTINUE
                        ISTART = ISTART + IADDT
  530                CONTINUE
                     ISTRTU = ISTRTU + IADDU
  520             CONTINUE
               END IF
            ELSE
               DO 550 IU = 2, JVAL, 2
                  UMIN1 = FLOAT(IU - 1)
                  ISTART = ISTRTU
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
                  DO 560 IT = 0, JVAL - IU, IPQ0X1
                     IPREV2 = ISTART + IADDPR
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                     DO 570 I = 1, ISTEPT
                        WORK1(ISTART + I) = UMIN1*WORK1(IPREV2 + I)
  570                CONTINUE
                     ISTART = ISTART + IADDT
  560             CONTINUE
                  ISTRTU = ISTRTU + IADDU
  550          CONTINUE
            END IF
C
C           R(T,U,V)
C
            ISTRTV = ISTRTJ + IADDV
            IADDPR = ISTEPJ - IADDV
            IF (PQZGT0) THEN
               ISTRTU = ISTRTV
               IUMAX = JVAL - 1
               DO 600 IU = 0, IUMAX, IPQ0Y1
                  ISTART = ISTRTU
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
                  DO 610 IT = 0, IUMAX - IU, IPQ0X1
                     IPREV1 = ISTART + IADDPR
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                     DO 620 I = 1, ISTEPT
                        WORK1(ISTART+I) = DPQZ(I)*WORK1(IPREV1+I)
  620                CONTINUE
                     ISTART = ISTART + IADDT
  610             CONTINUE
                  ISTRTU = ISTRTU + IADDU
  600          CONTINUE
               IF (JVAL .GT. 1) THEN
                  ISTRTV = ISTRTV + IADDV
                  DO 630 IV = 2,JVAL
                     ISTRTU = ISTRTV
                     VMIN1 = FLOAT(IV - 1)
                     IUMAX = JVAL - IV
                     DO 640 IU = 0, IUMAX, IPQ0Y1
                        ISTART = ISTRTU
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
                        DO 650 IT = 0, IUMAX - IU, IPQ0X1
                           IPREV1 = ISTART + IADDPR
                           IPREV2 = IPREV1 - IADDV
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                           DO 660 I = 1, ISTEPT
                              WORK1(ISTART+I) = DPQZ(I)*WORK1(IPREV1+I)
     *                                          + VMIN1*WORK1(IPREV2+I)
  660                      CONTINUE
                           ISTART = ISTART + IADDT
  650                   CONTINUE
                        ISTRTU = ISTRTU + IADDU
  640                CONTINUE
                     ISTRTV = ISTRTV + IADDV
  630             CONTINUE
               END IF
            ELSE
               DO 670 IV = 2, JVAL, 2
                  ISTRTU = ISTRTV
                  VMIN1 = FLOAT(IV - 1)
                  IUMAX = JVAL - IV
                  DO 680 IU = 0, IUMAX, IPQ0Y1
                     ISTART = ISTRTU
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
                     DO 690 IT = 0, IUMAX - IU, IPQ0X1
                        IPREV1 = ISTART + IADDPR
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                        DO 695 I = 1, ISTEPT
                           WORK1(ISTART+I) = VMIN1*WORK1(IPREV1+I)
  695                   CONTINUE
                        ISTART = ISTART + IADDT
  690                CONTINUE
                     ISTRTU = ISTRTU + IADDU
  680             CONTINUE
                  ISTRTV = ISTRTV + IADDV
  670          CONTINUE
            END IF
  200    CONTINUE
      END IF
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IPRINT = 0
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1010) JMAX
         WRITE (LUPRI, 1030) ISTEPT, ISTEPU, ISTEPV
         WRITE (LUPRI, 1040) NRTUV
         WRITE (LUPRI, 1050) PQXGT0, PQYGT0, PQZGT0
         WRITE (LUPRI, 1060) IWKHHS
         IF (IPRINT .GE. 20) THEN
            WRITE (LUPRI, 1100)
            DO 2000 IV = 0, JMAX, IPQ0Z + 1
               DO 2100 IU = 0, JMAX - IV, IPQ0Y + 1
                  DO 2200 IT = 0, JMAX - IU - IV, IPQ0X + 1
                     WRITE (LUPRI, 1110) IT, IU, IV
                     IADR = IWKHHS + IT*ISTEPT + IU*ISTEPU + IV*ISTEPV
                     WRITE (LUPRI, 1120) (WORK1(IADR+J),J=1,ISTEPT)
 2200             CONTINUE
 2100          CONTINUE
 2000       CONTINUE
         END IF
      END IF
      IPRINT = 0
      CALL HERDIS(ISTEPT,JMAX,NRTUV,IPRINT,WORK1,LWORK1)
      RETURN
C
 1000 FORMAT (//,' <<<<<<<<<< SUBROUTINE HERI >>>>>>>>>>')
 1010 FORMAT (/,'  JMAX    ',I7)
 1030 FORMAT (  '  ISTEP   ',3I7)
 1040 FORMAT (  '  NRTUV   ',I7)
 1050 FORMAT (  '  PQGT0   ',3L7)
 1060 FORMAT (  '  IWKHHS  ',I7)
 1100 FORMAT (//,' ***** HERMITIAN INTEGRALS ***** ')
 1110 FORMAT (//,'   R(', I2, ',', I2, ',', I2, ')')
 1120 FORMAT (/,/(1X,6F12.8))
C
C     ... end of HERI
C
      END
