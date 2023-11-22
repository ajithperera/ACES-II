      SUBROUTINE TWOODC(LA,LB,IOFF,NUCAB,ISTEPA,ISTEPB,SIGN,SIGNA,
     *                  SIGNB,IAB0,IADDAM,IADDA0,IADDAP,IADDBM,IADDB0,
     *                  IADDBP,DIFPA,DIFPB,HEXPPI,WORD,IPRINT)
C
C     TUH 84
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (ZERO =0.00 D00, ONE =1.00 D00, TWO =2.00 D00)
      CHARACTER WORD*4
      DIMENSION DIFPA(MXAOSQ), DIFPB(MXAOSQ), HEXPPI(MXAOSQ)
#include "lwork.par"
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
      ISTRT0 = IOFF
C
C     ****************************
C     ********** AB > 0 **********
C     ****************************
C
      IF (IAB0 .EQ. 0) THEN
C
C        ***** RUN OVER IA *****
C
         DO 100 IA = 0, LA
            ISTART = ISTRT0
C
C           ***** E(0,0) *****
C
            IF (IA .EQ. 0) THEN
               IADR = ISTART
               DO 200 I = 1, NUCAB
                  IADR = IADR + 1
                  WORK2(IADR) = ONE
  200          CONTINUE
C
C           ***** E(1,0) *****
C
            ELSE IF (IA .EQ. 1) THEN
               IADR = ISTART
               DO 300 I = 1, NUCAB
                  IADR = IADR + 1
                  WORK2(IADR) = DIFPA(I)
                  WORK2(IADR + NUCAB) = HEXPPI(I)
  300          CONTINUE
C
C           ***** E(2,0) *****
C
            ELSE IF (IA .EQ. 2) THEN
               IADR0 = ISTART
               IADR1 = IADR0 + NUCAB
               IADR2 = IADR1 + NUCAB
               DO 400 I = 1, NUCAB
                  IADR0 = IADR0 + 1
                  IADR1 = IADR1 + 1
                  IADR2 = IADR2 + 1
                  DIFPAI = DIFPA(I)
                  EXPPIH = HEXPPI(I)
                  WORK2(IADR0) = DIFPAI*DIFPAI + SIGN*EXPPIH
                  WORK2(IADR1) = TWO*DIFPAI*EXPPIH
                  WORK2(IADR2) = EXPPIH*EXPPIH
  400          CONTINUE
C
C            ***** E(IA,0) *****
C
            ELSE
               IADR0 = ISTART
               IADR2 = IADR0 + IA*NUCAB
               IADR1 = IADR2 - NUCAB
               DO 500 I = 1, NUCAB
                  IADR0 = IADR0 + 1
                  IADR1 = IADR1 + 1
                  IADR2 = IADR2 + 1
                  DIFPAI = DIFPA(I)
                  EXPPIH = HEXPPI(I)
                  WORK2(IADR0) = DIFPAI*WORK2(IADR0 + IADDA0)
     *                           + SIGN*WORK2(IADR0 + IADDAP)
                  WORK2(IADR1) = EXPPIH*WORK2(IADR1 + IADDAM)
     *                         + DIFPAI*WORK2(IADR1 + IADDA0)
                  WORK2(IADR2) = EXPPIH*WORK2(IADR2 + IADDAM)
  500          CONTINUE
               IADR = ISTART + NUCAB
               DO 510 IT = 1, IA - 2
                  T1 = SIGN*FLOAT(IT + 1)
                  DO 520 I = 1, NUCAB
                     IADR = IADR + 1
                     WORK2(IADR) = HEXPPI(I)*WORK2(IADR + IADDAM)
     *                           + DIFPA(I)*WORK2(IADR + IADDA0)
     *                           + T1*WORK2(IADR + IADDAP)
  520             CONTINUE
  510          CONTINUE
            END IF
            ISTART = ISTART + ISTEPB
C
C           ***** RUN OVER IB *****
C
            DO 600 IB = 1, LB
               IAB = IA + IB
C
C              ***** E(0,1) *****
C
               IF (IAB .EQ. 1) THEN
                  IADR = ISTART
                  DO 700 I = 1, NUCAB
                     IADR = IADR + 1
                     WORK2(IADR) = DIFPB(I)
                     WORK2(IADR + NUCAB) = HEXPPI(I)
  700             CONTINUE
               ELSE IF (IAB .EQ. 2) THEN
                  IADR0 = ISTART
                  IADR1 = IADR0 + NUCAB
                  IADR2 = IADR1 + NUCAB
C
C                 ***** E(0,2) *****
C
                  IF (IA .EQ. 0) THEN
                     DO 800 I = 1, NUCAB
                        IADR0 = IADR0 + 1
                        IADR1 = IADR1 + 1
                        IADR2 = IADR2 + 1
                        DIFPBI = DIFPB(I)
                        EXPPIH = HEXPPI(I)
                        WORK2(IADR0) = DIFPBI*DIFPBI + SIGN*EXPPIH
                        WORK2(IADR1) = TWO*DIFPBI*EXPPIH
                        WORK2(IADR2) = EXPPIH*EXPPIH
  800                CONTINUE
C
C                 ***** E(1,1) *****
C
                  ELSE
                     DO 810 I = 1, NUCAB
                        IADR0 = IADR0 + 1
                        IADR1 = IADR1 + 1
                        IADR2 = IADR2 + 1
                        DIFPAI = DIFPA(I)
                        DIFPBI = DIFPB(I)
                        EXPPIH = HEXPPI(I)
                        WORK2(IADR0) = DIFPAI*DIFPBI + SIGN*EXPPIH
                        WORK2(IADR1) = (DIFPAI + DIFPBI)*EXPPIH
                        WORK2(IADR2) = EXPPIH*EXPPIH
  810                CONTINUE
                  END IF
C
C              ***** E(IA,IB) *****
C
               ELSE
                  IADR0 = ISTART
                  IADR2 = IADR0 + IAB*NUCAB
                  IADR1 = IADR2 - NUCAB
                  DO 900 I = 1, NUCAB
                     IADR0 = IADR0 + 1
                     IADR1 = IADR1 + 1
                     IADR2 = IADR2 + 1
                     DIFPBI = DIFPB(I)
                     EXPPIH = HEXPPI(I)
                     WORK2(IADR0) = DIFPBI*WORK2(IADR0 + IADDB0)
     *                              + SIGN*WORK2(IADR0 + IADDBP)
                     WORK2(IADR1) = EXPPIH*WORK2(IADR1 + IADDBM)
     *                            + DIFPBI*WORK2(IADR1 + IADDB0)
                     WORK2(IADR2) = EXPPIH*WORK2(IADR2 + IADDBM)
  900             CONTINUE
                  IADR = ISTART + NUCAB
                  DO 910 IT = 1, IAB - 2
                     T1 = SIGN*FLOAT(IT + 1)
                     DO 920 I = 1, NUCAB
                        IADR = IADR + 1
                        WORK2(IADR) = HEXPPI(I)*WORK2(IADR + IADDBM)
     *                              + DIFPB(I)*WORK2(IADR + IADDB0)
     *                              + T1*WORK2(IADR + IADDBP)
  920                CONTINUE
  910             CONTINUE
               END IF
               ISTART = ISTART + ISTEPB
  600       CONTINUE
            ISTRT0 = ISTRT0 + ISTEPA
  100    CONTINUE
C
C     ****************************
C     ********** AB = 0 **********
C     ****************************
C
      ELSE
C
C        ***** RUN OVER IA *****
C
         DO 105 IA = 0, LA
            ISTART = ISTRT0
C
C           ***** E(0,0) *****
C
            IF (IA .EQ. 0) THEN
               IADR = ISTART
               DO 205 I = 1, NUCAB
                  IADR = IADR + 1
                  WORK2(IADR) = ONE
  205          CONTINUE
C
C           ***** E(1,0) *****
C
            ELSE IF (IA .EQ. 1) THEN
               IADR = ISTART + NUCAB
               DO 305 I = 1, NUCAB
                  IADR = IADR + 1
                  WORK2(IADR) = HEXPPI(I)
  305          CONTINUE
C
C           ***** E(2,0) *****
C
            ELSE IF (IA .EQ. 2) THEN
               IADR0 = ISTART
               IADR2 = IADR0 + NUCAB + NUCAB
               DO 405 I = 1, NUCAB
                  IADR0 = IADR0 + 1
                  IADR2 = IADR2 + 1
                  EXPPIH = HEXPPI(I)
                  WORK2(IADR0) = SIGN*EXPPIH
                  WORK2(IADR2) = EXPPIH*EXPPIH
  405          CONTINUE
C
C            ***** E(IA,0) *****
C
            ELSE
               IODDA = IBTAND(1,IA)
               IADR0 = ISTART
               IADR2 = IADR0 + IA*NUCAB
               IF (IODDA .EQ. 0) THEN
                  DO 505 I = 1, NUCAB
                     IADR0 = IADR0 + 1
                     IADR2 = IADR2 + 1
                     WORK2(IADR0) = SIGN*WORK2(IADR0 + IADDAP)
                     WORK2(IADR2) = HEXPPI(I)*WORK2(IADR2 + IADDAM)
  505             CONTINUE
               ELSE
                  DO 506 I = 1, NUCAB
                     IADR2 = IADR2 + 1
                     WORK2(IADR2) = HEXPPI(I)*WORK2(IADR2 + IADDAM)
  506             CONTINUE
               END IF
               DO 515 IT = 2 - IODDA, IA - 2, 2
                  IADR = ISTART + IT*NUCAB
                  T1 = SIGN*FLOAT(IT + 1)
                  DO 525 I = 1, NUCAB
                     IADR = IADR + 1
                     WORK2(IADR) = HEXPPI(I)*WORK2(IADR + IADDAM)
     *                           + T1*WORK2(IADR + IADDAP)
  525             CONTINUE
  515          CONTINUE
            END IF
            ISTART = ISTART + ISTEPB
C
C           ***** RUN OVER IB *****
C
            DO 605 IB = 1, LB
               IAB = IA + IB
C
C              ***** E(0,1) *****
C
               IF (IAB .EQ. 1) THEN
                  IADR = ISTART + NUCAB
                  DO 705 I = 1, NUCAB
                     IADR = IADR + 1
                     WORK2(IADR) = HEXPPI(I)
  705             CONTINUE
C
C              ***** E(1,1) AND E(0,2) *****
C
               ELSE IF (IAB .EQ. 2) THEN
                  IADR0 = ISTART
                  IADR2 = IADR0 + NUCAB + NUCAB
                  DO 805 I = 1, NUCAB
                     IADR0 = IADR0 + 1
                     IADR2 = IADR2 + 1
                     EXPPIH = HEXPPI(I)
                     WORK2(IADR0) = SIGN*EXPPIH
                     WORK2(IADR2) = EXPPIH*EXPPIH
  805             CONTINUE
C
C              ***** E(IA,IB) *****
C
               ELSE
                  IODDAB = IBTAND(1,IAB)
                  IADR0 = ISTART
                  IADR2 = IADR0 + IAB*NUCAB
                  IF (IODDAB .EQ. 0) THEN
                     DO 905 I = 1, NUCAB
                        IADR0 = IADR0 + 1
                        IADR2 = IADR2 + 1
                        WORK2(IADR0) = SIGN*WORK2(IADR0 + IADDBP)
                        WORK2(IADR2) = HEXPPI(I)*WORK2(IADR2 + IADDBM)
  905                CONTINUE
                  ELSE
                     DO 906 I = 1, NUCAB
                        IADR2 = IADR2 + 1
                        WORK2(IADR2) = HEXPPI(I)*WORK2(IADR2 + IADDBM)
  906                CONTINUE
                  END IF
                  DO 915 IT = 2 - IODDAB, IAB - 2, 2
                     IADR = ISTART + IT*NUCAB
                     T1 = SIGN*FLOAT(IT + 1)
                     DO 925 I = 1, NUCAB
                        IADR = IADR + 1
                        WORK2(IADR) = HEXPPI(I)*WORK2(IADR + IADDBM)
     *                              + T1*WORK2(IADR + IADDBP)
  925                CONTINUE
  915             CONTINUE
               END IF
               ISTART = ISTART + ISTEPB
  605       CONTINUE
            ISTRT0 = ISTRT0 + ISTEPA
  105    CONTINUE
      END IF
      IF (IPRINT .LT. 10) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      WRITE (LUPRI, 1000)
      WRITE (LUPRI, 1010) LA, LB
      WRITE (LUPRI, 1020) IOFF
      WRITE (LUPRI, 1030) IAB0
      IF (IPRINT .LT. 20) RETURN
      IADR0 = IOFF
      DO 2000 IA = 0, LA
         IADR = IADR0
         DO 2100 IB = 0, LB
            IADRT = IADR
            DO 2200 IT = 0, IA + IB
               IODD = IBTAND(IA + IB - IT,IAB0)
               IF (IODD .EQ. 0) THEN
                  WRITE (LUPRI, 1100) WORD, IA, IB, IT
                  WRITE (LUPRI, 1130) (WORK2(IADRT + I), I = 1, NUCAB)
               END IF
               IADRT = IADRT + NUCAB
 2200       CONTINUE
            IADR = IADR + ISTEPB
 2100    CONTINUE
         IADR0 = IADR0 + ISTEPA
 2000 CONTINUE
      RETURN
 1000 FORMAT (/,'  <<<<<<<<<< SUBROUTINE TWOODC >>>>>>>>>>',/)
 1010 FORMAT ('  LA/B:     ',2I7)
 1020 FORMAT ('  IOFF:     ',I7)
 1030 FORMAT ('  IAB0:     ',I7)
 1100 FORMAT (/,1X,A4,'(',I1,',',I1,';',I1,')',/)
 1130 FORMAT(1X,6F12.8)
      END
