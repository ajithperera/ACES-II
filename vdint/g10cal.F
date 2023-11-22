      SUBROUTINE G10CAL(WORK1,LWORK1)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (AZERO = 0.00 D00)
C HCGT0 should be (2*mxang+1)*(2*mxang+2)*(2*mxang+3)/6, but since we
C don't have a parameter for mxang (=5) yet, let's hard code it.
      LOGICAL HCGT0(1000), PATH1, SEGCON, TPRI12, TCON12
      DIMENSION WORK1(1)
#include "lwork.par"
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /CWORK3/ WK3LOW, WORK3(LWORK3), WK3HGH
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      COMMON /PRIVEC/ NPRM12(MXAOSQ), NPRM34(MXAOSQ),
     *                LSTP12(MXAOSQ), LSTP34(MXAOSQ),
     *                LSTO12(MXAOSQ), LSTO34(MXAOSQ),
     *                LASTP(MXAOSQ),  LASTO(MXAOSQ)
#include "ninahh.par"
      COMMON /CINAHH/ INAHH(NINAHH)
      COMMON /CC1INF/ LVAL12, MVAL12, NVAL12, INCRMT, INCRMU, INCRMV,
     *                ISTRET, ISTREU, ISTREV, ISTEPT, ISTEPU, ISTEPV,
     *                IPQ0X, IPQ0Y, IPQ0Z, NUC1, NUC2, NUC12, NUC34,
     *                NCCPP, MAX34, NORB1, NORB2, NORB12,
     *                TPRI12, TCON12,KHKT12,
     *                INTHC0, IPRINT, PATH1, JSTRH(969),
     *                JODDH(969), IADDPV, ISTRCF, IOFFVC, ICMP12,
     *                MAXDER, NSET1, NSET2,
     *                NWORK0, NWORKH, NWORKE, SEGCON
      COMMON /DHCINF/ IDHC(10), NDHC(10)
      COMMON /ODCADR/ IE1X00, IE1X10, IE1X01, IE1X20, IE1X11, IE1X02,
     *                IE1Y00, IE1Y10, IE1Y01, IE1Y20, IE1Y11, IE1Y02,
     *                IE1Z00, IE1Z10, IE1Z01, IE1Z20, IE1Z11, IE1Z02,
     *                IE2X00, IE2X10, IE2X01, IE2X20, IE2X11, IE2X02,
     *                IE2Y00, IE2Y10, IE2Y01, IE2Y20, IE2Y11, IE2Y02,
     *                IE2Z00, IE2Z10, IE2Z01, IE2Z20, IE2Z11, IE2Z02,
     *                IELAST, ILST12
      COMMON /INT10/ NR34, NRSTEP, IOFFHC, IEXADR, IEYADR, IEZADR
      SAVE INTHCX,INTHCY,INTHCZ
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C  EVALUATION OF INTEGRALS 
C
C  CARTESIAN INTEGRALS 
C
C  SPECIAL CASE: (SS|XY) 
C
      IF (LVAL12 + MVAL12 + NVAL12 .EQ. 0) THEN
       NC134 = NUC1*NUC34
       IF (TPRI12) THEN
        LSCR1 = NUC2*NC134
       ELSE
        LSCR1 = NORB2*NC134
       END IF
       LSCR2 = NORB2*NC134
       IF (LSCR1 + LSCR2 .GT. LWORK3) THEN
        LSCR12 = LSCR1 + LSCR2
        WRITE (LUPRI,'(/A,/,A,I6,A,I6,A)')
     *       ' Work space for CNTOUT insufficient.',
     *       ' Current value of LWORK3',LWORK3,
     *       ' must be increased to at least ',LSCR12
        CALL ERREX
       END IF
       ISCR1 = 1
       ISCR2 = ISCR1 + LSCR2
       DO 100 ITUV = 1, NR34
        IF (JODDH(ITUV) .EQ. 0) THEN
         INTHH  = INAHH(JSTRH(ITUV)) + 1
         INTHCX = IOFFHC + (ITUV - 1)*NCCPP + 1
         CALL CNTOUT(WORK1(INTHH),WORK1(INTHCX),
     *               WORK3(ISCR1),WORK3(ISCR2),LSCR1,LSCR2,
     *               ISTRCF,IOFFVC,
     *               NSET2,NSET1,NUC2,NUC1,NORB2,NORB1,
     *               NUC34,TPRI12,TCON12)
        END IF
  100  CONTINUE
C
C  GENERAL CASE: (XY|ZW) 
C
      ELSE
C
C  FIRST CALCULATE PRIMITIVE INTEGRALS
C
       IOFFEX = IEXADR + ISTRET
       IOFFEY = IEYADR + ISTREU
       IOFFEZ = IEZADR + ISTREV
       MAXT = LVAL12
       MAXU = MVAL12
       MAXV = NVAL12
       MINT = IBTAND(MAXT,INCRMT - 1)
       MINU = IBTAND(MAXU,INCRMU - 1)
       MINV = IBTAND(MAXV,INCRMV - 1)
       DO 200 I = 1, NR34
        HCGT0(I) = .FALSE.
  200  CONTINUE
       DO 210 I = 1, NWORK0
        WORK1(IWKLST + I) = AZERO
  210  CONTINUE
       DO 300 IV = MINV, MAXV, INCRMV
        JSTREV = IOFFEZ + IV*NUC12
        ISTRHV = IV*ISTEPV
        IODDZ  = IBTSHL(IBTAND(IPQ0Z,IV),2)
        DO 310 IU = MINU, MAXU, INCRMU
         JSTREU = IOFFEY + IU*NUC12
         ISTRHU = ISTRHV + IU*ISTEPU
         IODDYZ = IBTOR(IODDZ,IBTSHL(IBTAND(IPQ0Y,IU),1))
         DO 320 IT = MINT, MAXT, INCRMT
          JSTRET = IOFFEX + IT*NUC12
          IODXYZ = IBTOR(IODDYZ,IBTAND(IPQ0X,IT))
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
          DO 400 I = 1, NUC12
           ECOEFI = WORK2(JSTRET + I)
     *            * WORK2(JSTREU + I)
     *            * WORK2(JSTREV + I)
           IJ = (I - 1)*NUC34
           DO 410 J = 1, NUC34
            WORK3(IJ + J) = ECOEFI
  410      CONTINUE
  400     CONTINUE
          INTHH0 = ISTRHU + IT
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
          DO 500 ITUV = 1, NR34
           IF (JODDH(ITUV) .EQ. IODXYZ) THEN
            HCGT0(ITUV) = .TRUE.
            IWORK = IWKLST + (ITUV - 1)*ISTEPT
            INTHH = INAHH(INTHH0 + JSTRH(ITUV))
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
*VOCL LOOP,NOVREC
CDIR$ IVDEP
            DO 510 I = 1, ISTEPT
             WORK1(IWORK + I) = WORK1(IWORK + I)
     *              + WORK3(I)*WORK1(INTHH + I)
  510       CONTINUE
           END IF
  500     CONTINUE
  320    CONTINUE
  310   CONTINUE
  300  CONTINUE
C
C  NOW CALCULATE CONTRACTED INTEGRALS
C
       NC134 = NUC1*NUC34
       IF (TPRI12) THEN
        LSCR1 = NUC2*NC134
       ELSE
        LSCR1 = NORB2*NC134
       END IF
       LSCR2 = NORB2*NC134
       IF (LSCR1 + LSCR2 .GT. LWORK3) THEN
        LSCR12 = LSCR1 + LSCR2
        WRITE (LUPRI,'(/A,/,A,I6,A,I6,A)')
     *           ' Work space for CNTOUT insufficient.',
     *           ' Current value of LWORK3',LWORK3,
     *           ' must be increased to at least ',LSCR12
         CALL ERREX
        END IF
        ISCR1 = 1
        ISCR2 = ISCR1 + LSCR2
        DO 600 ITUV = 1, NR34
         IF (HCGT0(ITUV)) THEN
          IWORK0 = IWKLST + (ITUV - 1)*ISTEPT + 1
          INTHCX = IOFFHC + (ITUV - 1)*NCCPP  + 1
          CALL CNTOUT(WORK1(IWORK0),WORK1(INTHCX),
     *                WORK3(ISCR1),WORK3(ISCR2),LSCR1,LSCR2,
     *                ISTRCF,IOFFVC,
     *                NSET2,NSET1,NUC2,NUC1,NORB2,NORB1,
     *                NUC34,TPRI12,TCON12)
         END IF
  600   CONTINUE
      END IF
C
      IOFFHC = IOFFHC + NR34*NCCPP
      RETURN
 2000 FORMAT (//,1X,'<<<<<<<<<< SUBROUTINE G10CAL >>>>>>>>>>',/)
 2010 FORMAT (1X,'IOFFEX/Y/Z/ ',3I7)
 2020 FORMAT (1X,'LOOP PARAMETERS X-DIR ',3I7)
 2030 FORMAT (1X,'LOOP PARAMETERS Y-DIR ',3I7)
 2040 FORMAT (1X,'LOOP PARAMETERS Z-DIR ',3I7)
 2050 FORMAT (/,1X,'***** HERMITIAN-CARTESIAN INTEGRALS *****',/)
 3000 FORMAT (1X,'ITUV   ',I7)
 3030 FORMAT (1X,'ORB12  ',I7/,(1P,6E12.4))
      END
