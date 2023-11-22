      SUBROUTINE C1HCAL(WORK1,LWORK1)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (ZERO = 0.00 D00, ONE = 1.00 D00)
C HCGT0 should be (2*mxang+1)*(2*mxang+2)*(2*mxang+3)/6, but since we
C don't have a parameter for mxang (=5) yet, let's hard code it.
      LOGICAL HCGT0(1000), PATH1, SEGCON, TPRI12,TCON12,
     *        DHCHX, DHCHY, DHCHZ, DHCEX, DHCEY, DHCEZ,
     *        DHCEX1, DHCEX2, DHCEY1, DHCEY2, DHCEZ1, DHCEZ2
      DIMENSION NR34H(4), IADDHH(4), JODDIF(4), IOFFHC(4),
     *          SGN(4)
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
     *                TPRI12, TCON12, KHKT12,
     *                INTHC0, IPRINT, PATH1, JSTRH(969),
     *                JODDH(969), IADDPV, ISTRCF, IOFFVC, ICMP12,
     *                MAXDER, NSET1, NSET2,
     *                NWORK0, NWORKH, NWORKE, SEGCON
      COMMON /ODCADR/ IE1X00, IE1X10, IE1X01, IE1X20, IE1X11, IE1X02,
     *                IE1Y00, IE1Y10, IE1Y01, IE1Y20, IE1Y11, IE1Y02,
     *                IE1Z00, IE1Z10, IE1Z01, IE1Z20, IE1Z11, IE1Z02,
     *                IE2X00, IE2X10, IE2X01, IE2X20, IE2X11, IE2X02,
     *                IE2Y00, IE2Y10, IE2Y01, IE2Y20, IE2Y11, IE2Y02,
     *                IE2Z00, IE2Z10, IE2Z01, IE2Z20, IE2Z11, IE2Z02,
     *                IELAST, ILST12
      COMMON /CRSDIR/ DHCHX, DHCHY, DHCHZ,
     *                DHCEX, DHCEX1, DHCEX2,
     *                DHCEY, DHCEY1, DHCEY2,
     *                DHCEZ, DHCEZ1, DHCEZ2
      COMMON /DHCINF/ IHC00,
     *                IHCHX, IHCEX1, IHCEX2,
     *                IHCHY, IHCEY1, IHCEY2,
     *                IHCHZ, IHCEZ1, IHCEZ2,
     *                NHC00,
     *                NHCHX, NHCEX1, NHCEX2,
     *                NHCHY, NHCEY1, NHCEY2,
     *                NHCHZ, NHCEZ1, NHCEZ2
      COMMON/INT1H/SGN,NR34H,IADDHH,JODDIF,IOFFHC,NHCGT0,
     *             IEXADR,IEYADR,IEZADR,NCOOR
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     ***** CALCULATION OF INTEGRALS - ENTRY POINT C1HCAL *****
C
C
C     ***** PRIMITIVE CARTESIAN INTEGRALS *****
C
      IPRINT = 21
      IOFFEX = IEXADR + ISTRET
      IOFFEY = IEYADR + ISTREU
      IOFFEZ = IEZADR + ISTREV
      MAXT = LVAL12
      MAXU = MVAL12
      MAXV = NVAL12
      MINT = IBTAND(MAXT,INCRMT - 1)
      MINU = IBTAND(MAXU,INCRMU - 1)
      MINV = IBTAND(MAXV,INCRMV - 1)
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 2000)
         WRITE (LUPRI, 2010) IOFFEX, IOFFEY, IOFFEZ
         WRITE (LUPRI, 2020) MINT, MAXT, INCRMT
         WRITE (LUPRI, 2030) MINU, MAXU, INCRMU
         WRITE (LUPRI, 2040) MINV, MAXV, INCRMV
         IF (IPRINT .GE. 20) WRITE (LUPRI, 2050)
      END IF
      DO 100 I = 1, NHCGT0
         HCGT0(I) = .FALSE.
  100 CONTINUE
      DO 150 I = 1, NWORKH
         WORK1(IWKLST + I) = ZERO
  150 CONTINUE
      DO 200 IV = MINV, MAXV, INCRMV
         JSTREV = IOFFEZ + IV*NUC12
         ISTRHV = IV*ISTEPV
         IODDZ = IBTSHL(IBTAND(IPQ0Z,IV),2)
         DO 210 IU = MINU, MAXU, INCRMU
            JSTREU = IOFFEY + IU*NUC12
            ISTRHU = ISTRHV + IU*ISTEPU
            IODDYZ = IBTOR(IODDZ,IBTSHL(IBTAND(IPQ0Y,IU),1))
            DO 220 IT = MINT, MAXT, INCRMT
               JSTRET = IOFFEX + IT*NUC12
               ISTRTH = ISTRHU + IT
               IODXYZ = IBTOR(IODDYZ,IBTAND(IPQ0X,IT))
               IJ = 0
               DO 300 I = 1, NUC12
                  ECOEFI = WORK2(JSTRET + I)
     *                   * WORK2(JSTREU + I)
     *                   * WORK2(JSTREV + I)
                  DO 310 J = 1, NUC34
                     IJ = IJ + 1
                     WORK3(IJ) = ECOEFI
  310             CONTINUE
  300          CONTINUE
               ITUV0 = 0
               IWORK = IWKLST
               DO 400 ICOOR = 1, NCOOR
                  NR34 = NR34H(ICOOR)
                  INTHH0 = ISTRTH + IADDHH(ICOOR)
                  JODXYZ = IBTXOR(IODXYZ,JODDIF(ICOOR))
                  DO 500 ITUV = 1, NR34
                     IF (JODDH(ITUV) .EQ. JODXYZ) THEN
                        HCGT0(ITUV0 + ITUV) = .TRUE.
                        INTHH = INAHH(INTHH0 + JSTRH(ITUV))
                        DO 510 I = 1, ISTEPT
                           WORK1(IWORK + I) = WORK1(IWORK + I)
     *                             + WORK3(I)*WORK1(INTHH + I)
  510                   CONTINUE
                     END IF
                     IWORK = IWORK + ISTEPT
  500             CONTINUE
                  ITUV0 = ITUV0 + NR34
  400          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
C     ******************************************
C     ***** CONTRACTED CARTESIAN INTEGRALS *****
C     ******************************************
C
      IWORK0 = IWKLST
      ITUV0 = 0
      DO 600 ICOOR = 1, NCOOR
         INTHC0 = IOFFHC(ICOOR)
         SIGNI = SGN(ICOOR)
         NR34 = NR34H(ICOOR)
         DO 700 ITUV = 1, NR34
            IF (HCGT0(ITUV0 + ITUV)) THEN
               INTHC = INTHC0
               IWORK1 = IWORK0
               DO 800 I = 1, NORB12
                  KPRM12 = NPRM12(I + IADDPV)
                  DO 810 J = 1, NUC34
                     SUM = ZERO
                     IWORK = IWORK1 + J
                     DO 820 K = 1, KPRM12
                        SUM = SUM + WORK1(IWORK)
                        IWORK = IWORK + NUC34
  820                CONTINUE
                     INTHC = INTHC + 1
                     WORK1(INTHC) = SIGNI*SUM
  810             CONTINUE
                  IWORK1 = IWORK1 + NUC34*KPRM12
  800          CONTINUE
               IF (IPRINT .GE. 20) THEN
                  WRITE (LUPRI, 3000) ICOOR
                  WRITE (LUPRI, 3010) ITUV
                  WRITE (LUPRI, 3020) IWORK0
                  WRITE (LUPRI, 3030) INTHC0
                  INTHC = INTHC0
                  DO 900 I = 1, NORB12
                     WRITE (LUPRI, 3040) I
                     WRITE (LUPRI, 3050)
     *                     (WORK1(INTHC + J), J = 1, NUC34)
                     INTHC = INTHC + NUC34
  900             CONTINUE
               END IF
            END IF
            INTHC0 = INTHC0 + NCCPP
            IWORK0 = IWORK0 + ISTEPT
  700    CONTINUE
         ITUV0 = ITUV0 + NR34
         IOFFHC(ICOOR) = INTHC0
  600 CONTINUE
      RETURN
 2000 FORMAT (//,1X,'<<<<<<<<<< SUBROUTINE C1HINT >>>>>>>>>>')
 2010 FORMAT (1X,'IOFFEX/Y/Z/ ',3I7)
 2020 FORMAT (1X,'LOOP PARAMETERS X-DIR ',3I7)
 2030 FORMAT (1X,'LOOP PARAMETERS Y-DIR ',3I7)
 2040 FORMAT (1X,'LOOP PARAMETERS Z-DIR ',3I7)
 2050 FORMAT (//,1X,'***** HERMITIAN-CARTESIAN INTEGRALS *****',/)
 3000 FORMAT (1X,'ICOOR  ',I7)
 3010 FORMAT (1X,'ITUV   ',I7)
 3020 FORMAT (1X,'IWORK0 ',I7)
 3030 FORMAT (1X,'INTHC0 ',I7)
 3040 FORMAT (1X,'ORB12  ',I7)
 3050 FORMAT (1P,6E12.4)
      END
