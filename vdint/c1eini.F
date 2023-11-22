      SUBROUTINE C1EINI
C
C     TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (AZERO = 0.00 D00)
C HCGT0 should be (2*mxang+1)*(2*mxang+2)*(2*mxang+3)/6, but since we
C don't have a parameter for mxang (=5) yet, let's hard code it.
      LOGICAL HCGT0(1000), PATH1, SEGCON, TPRI12,TCON12,
     *        DHCHX, DHCHY, DHCHZ, DHCEX, DHCEY, DHCEZ,
     *        DHCEX1, DHCEX2, DHCEY1, DHCEY2, DHCEZ1, DHCEZ2
      DIMENSION JOFFEX(6), JOFFEY(6), JOFFEZ(6), INCMXT(3),
     *          INCMXU(3), INCMXV(3), JSTRAT(3), JENDAT(3),
     *          IOFFHC(6)
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      COMMON /PRIVEC/ NPRM12(MXAOSQ), NPRM34(MXAOSQ),
     *                LSTP12(MXAOSQ), LSTP34(MXAOSQ),
     *                LSTO12(MXAOSQ), LSTO34(MXAOSQ),
     *                LASTP(MXAOSQ),  LASTO(MXAOSQ)
#include "ninahh.par"
      COMMON /CINAHH/ INAHH(NINAHH)
      COMMON /CCFCOM/ CONT1 (MXCONT*MXAOVC), CONT2 (MXCONT*MXAOVC),
     *                CONT3 (MXCONT*MXAOVC), CONT4 (MXCONT*MXAOVC),
     *                CONTT1(MXCONT*MXAOVC), CONTT2(MXCONT*MXAOVC),
     *                CONTT3(MXCONT*MXAOVC), CONTT4(MXCONT*MXAOVC),
     *                NUC1X (MXAOVC),        NUC2X (MXAOVC),
     *                NUC3X (MXAOVC),        NUC4X (MXAOVC),
     *                NRC1X (MXAOVC),        NRC2X (MXAOVC),
     *                NRC3X (MXAOVC),        NRC4X (MXAOVC)
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
      COMMON /DHCINF/ IHC00,
     *                IHCHX, IHCEX1, IHCEX2,
     *                IHCHY, IHCEY1, IHCEY2,
     *                IHCHZ, IHCEZ1, IHCEZ2,
     *                NHC00,
     *                NHCHX, NHCEX1, NHCEX2,
     *                NHCHY, NHCEY1, NHCEY2,
     *                NHCHZ, NHCEZ1, NHCEZ2
      COMMON /CRSDIR/ DHCHX, DHCHY, DHCHZ,
     *                DHCEX, DHCEX1, DHCEX2,
     *                DHCEY, DHCEY1, DHCEY2,
     *                DHCEZ, DHCEZ1, DHCEZ2
      COMMON/INT1E/NCOOR,NR34,JOFFEX,JOFFEY,JOFFEZ,INCMXT,INCMXU,
     &             INCMXV,JSTRAT, JENDAT, IOFFHC,NRCCPP
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     ***********************************************
C     ***** INITIALIZATION - ENTRY POINT C1EINI *****
C     ***********************************************
C
      IF (PATH1) THEN
         IEX00 = IE1X00
         IEX10 = IE1X10
         IEX01 = IE1X01
         IEY00 = IE1Y00
         IEY10 = IE1Y10
         IEY01 = IE1Y01
         IEZ00 = IE1Z00
         IEZ10 = IE1Z10
         IEZ01 = IE1Z01
      ELSE
         IEX00 = IE2X00
         IEX10 = IE2X10
         IEX01 = IE2X01
         IEY00 = IE2Y00
         IEY10 = IE2Y10
         IEY01 = IE2Y01
         IEZ00 = IE2Z00
         IEZ10 = IE2Z10
         IEZ01 = IE2Z01
      END IF
      ICOOR = 0
      IATOM = 0
      INTHC = IWKLST
      MAX34E = MAX34 + MAXDER - 1
      NR34 = (MAX34E + 1)*(MAX34E + 2)*(MAX34E + 3)/6
      NRCCPP = NR34*NCCPP
      KHNRCP = KHKT12*NRCCPP
C
C     ***** X-DIRECTION *****
C
      IF (DHCEX) THEN
         ICOOR = ICOOR + 1
         JSTRAT(ICOOR) = IATOM + 1
         INCMXT(ICOOR) = 1
         INCMXU(ICOOR) = 0
         INCMXV(ICOOR) = 0
         IF (DHCEX1) THEN
            IATOM = IATOM + 1
            JOFFEX(IATOM) = IEX10
            JOFFEY(IATOM) = IEY00
            JOFFEZ(IATOM) = IEZ00
            IOFFHC(IATOM) = INTHC
            IHCEX1 = INTHC
            NHCEX1 = NRCCPP
            INTHC = INTHC + KHNRCP
         END IF
         IF (DHCEX2) THEN
            IATOM = IATOM + 1
            JOFFEX(IATOM) = IEX01
            JOFFEY(IATOM) = IEY00
            JOFFEZ(IATOM) = IEZ00
            IOFFHC(IATOM) = INTHC
            IHCEX2 = INTHC
            NHCEX2 = NRCCPP
            INTHC = INTHC + KHNRCP
         END IF
         JENDAT(ICOOR) = IATOM
      END IF
C
C     ***** Y-DIRECTION *****
C
      IF (DHCEY) THEN
         ICOOR = ICOOR + 1
         JSTRAT(ICOOR) = IATOM + 1
         INCMXT(ICOOR) = 0
         INCMXU(ICOOR) = 1
         INCMXV(ICOOR) = 0
         IF (DHCEY1) THEN
            IATOM = IATOM + 1
            JOFFEX(IATOM) = IEX00
            JOFFEY(IATOM) = IEY10
            JOFFEZ(IATOM) = IEZ00
            IOFFHC(IATOM) = INTHC
            IHCEY1 = INTHC
            NHCEY1 = NRCCPP
            INTHC = INTHC + KHNRCP
         END IF
         IF (DHCEY2) THEN
            IATOM = IATOM + 1
            JOFFEX(IATOM) = IEX00
            JOFFEY(IATOM) = IEY01
            JOFFEZ(IATOM) = IEZ00
            IOFFHC(IATOM) = INTHC
            IHCEY2 = INTHC
            NHCEY2 = NRCCPP
            INTHC = INTHC + KHNRCP
         END IF
         JENDAT(ICOOR) = IATOM
      END IF
C
C     ***** Z-DIRECTION *****
C
      IF (DHCEZ) THEN
         ICOOR = ICOOR + 1
         JSTRAT(ICOOR) = IATOM + 1
         INCMXT(ICOOR) = 0
         INCMXU(ICOOR) = 0
         INCMXV(ICOOR) = 1
         IF (DHCEZ1) THEN
            IATOM = IATOM + 1
            JOFFEX(IATOM) = IEX00
            JOFFEY(IATOM) = IEY00
            JOFFEZ(IATOM) = IEZ10
            IOFFHC(IATOM) = INTHC
            IHCEZ1 = INTHC
            NHCEZ1 = NRCCPP
            INTHC = INTHC + KHNRCP
         END IF
         IF (DHCEZ2) THEN
            IATOM = IATOM + 1
            JOFFEX(IATOM) = IEX00
            JOFFEY(IATOM) = IEY00
            JOFFEZ(IATOM) = IEZ01
            IOFFHC(IATOM) = INTHC
            IHCEZ2 = INTHC
            NHCEZ2 = NRCCPP
            INTHC = INTHC + KHNRCP
         END IF
         JENDAT(ICOOR) = IATOM
      END IF
      IWKLST = INTHC
      NCOOR = ICOOR
      NWORKE = NR34*ISTEPT
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1020) DHCEX, DHCEX1, DHCEX2
         WRITE (LUPRI, 1030) DHCEY, DHCEY1, DHCEY2
         WRITE (LUPRI, 1040) DHCEZ, DHCEZ1, DHCEZ2
         WRITE (LUPRI, 1045) NWORKE
         WRITE (LUPRI, 1050) NCOOR
         WRITE (LUPRI, 1060) IATOM
         WRITE (LUPRI, 1070) NR34
         WRITE (LUPRI, 1080) (INCMXT(I), I = 1, NCOOR)
         WRITE (LUPRI, 1090) (INCMXU(I), I = 1, NCOOR)
         WRITE (LUPRI, 1100) (INCMXV(I), I = 1, NCOOR)
         WRITE (LUPRI, 1120) (JOFFEX(I), I = 1, IATOM)
         WRITE (LUPRI, 1130) (JOFFEY(I), I = 1, IATOM)
         WRITE (LUPRI, 1140) (JOFFEZ(I), I = 1, IATOM)
      END IF
      RETURN
 1000 FORMAT (//,1X,'***** SUBROUTINE C1EINT - INITIALIZATION *****',/)
 1020 FORMAT (1X,'DHCEX(1/2)',3L7)
 1030 FORMAT (1X,'DHCEY(1/2)',3L7)
 1040 FORMAT (1X,'DHCEZ(1/2)',3L7)
 1045 FORMAT (1X,'NWORKE    ',I7)
 1050 FORMAT (1X,'NCOOR     ',I7)
 1060 FORMAT (1X,'IATOM     ',I7)
 1070 FORMAT (1X,'NR34      ',I7)
 1080 FORMAT (1X,'INCMXT    ',(4I7))
 1090 FORMAT (1X,'INCMXU    ',(4I7))
 1100 FORMAT (1X,'INCMXV    ',(4I7))
 1120 FORMAT (1X,'JOFFEX    ',(7I7))
 1130 FORMAT (1X,'JOFFEY    ',(7I7))
 1140 FORMAT (1X,'JOFFEZ    ',(7I7))
      END
