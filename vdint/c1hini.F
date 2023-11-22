      SUBROUTINE C1HINI
C
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
     &             IEXADR,IEYADR,IEZADR,NCOOR
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     ***********************************************
C     ***** INITIALIZATION - ENTRY POINT C1HINI *****
C     ***********************************************
C
      IF (PATH1) THEN
         IEXADR = IE1X00
         IEYADR = IE1Y00
         IEZADR = IE1Z00
         SIGNI = ONE
      ELSE
         IEXADR = IE2X00
         IEYADR = IE2Y00
         IEZADR = IE2Z00
         SIGNI = - ONE
      END IF
      ICOOR = 0
      INTHC = IWKLST
      NRSUM = 0
C
C     ***** UNDIFFERENTIATED INTEGRALS *****
C
      ICOOR = ICOOR + 1
      IADDHH(ICOOR) = 0
      JODDIF(ICOOR) = 0
      MAX34H = MAX34 + MAXDER
      NR34 = (MAX34H + 1)*(MAX34H + 2)*(MAX34H + 3)/6
      NRCCPP = NR34*NCCPP
      NR34H(ICOOR) = NR34
      SGN(ICOOR) = ONE
      IOFFHC(ICOOR) = INTHC
      IHC00 = INTHC
      NHC00 = NRCCPP
      INTHC = INTHC + KHKT12*NRCCPP
      NRSUM = NRSUM + NR34
C
C     ***** DIFFERENTIATED INTEGRALS *****
C
      MAX34H = MAX34 + MAXDER - 1
      NR34 = (MAX34H + 1)*(MAX34H + 2)*(MAX34H + 3)/6
      NRCCPP = NR34*NCCPP
      KHNRCP = KHKT12*NRCCPP
C
C     ***** X-DIRECTION *****
C
      IF (DHCHX) THEN
         ICOOR = ICOOR + 1
         IADDHH(ICOOR) = 1
         JODDIF(ICOOR) = IPQ0X
         NR34H(ICOOR) = NR34
         SGN(ICOOR) = SIGNI
         IOFFHC(ICOOR) = INTHC
         IHCHX = INTHC
         NHCHX = NRCCPP
         INTHC = INTHC + KHNRCP
         NRSUM = NRSUM + NR34
      END IF
C
C     ***** Y-DIRECTION *****
C
      IF (DHCHY) THEN
         ICOOR = ICOOR + 1
         IADDHH(ICOOR) = ISTEPU
         JODDIF(ICOOR) = IBTSHL(IPQ0Y,1)
         NR34H(ICOOR) = NR34
         SGN(ICOOR) = SIGNI
         IOFFHC(ICOOR) = INTHC
         IHCHY = INTHC
         NHCHY = NRCCPP
         INTHC = INTHC + KHNRCP
         NRSUM = NRSUM + NR34
      END IF
C
C     ***** Z-DIRECTION *****
C
      IF (DHCHZ) THEN
         ICOOR = ICOOR + 1
         IADDHH(ICOOR) = ISTEPV
         JODDIF(ICOOR) = IBTSHL(IPQ0Z,2)
         NR34H(ICOOR) = NR34
         SGN(ICOOR) = SIGNI
         IOFFHC(ICOOR) = INTHC
         IHCHZ = INTHC
         NHCHZ = NRCCPP
         INTHC = INTHC + KHNRCP
         NRSUM = NRSUM + NR34
      END IF
      NCOOR = ICOOR
      IWKLST = INTHC
      NHCGT0 = NRSUM
      NWORKH = NRSUM*ISTEPT
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 1000)
         WRITE (LUPRI, 1010) IEXADR, IEYADR, IEZADR
         WRITE (LUPRI, 1020) NCOOR
         WRITE (LUPRI, 1030) NWORKH
         WRITE (LUPRI, 1040) NHCGT0
         WRITE (LUPRI, 1050) (IADDHH(I), I = 1, NCOOR)
         WRITE (LUPRI, 1060) (JODDIF(I), I = 1, NCOOR)
         WRITE (LUPRI, 1070) (IOFFHC(I), I = 1, NCOOR)
         WRITE (LUPRI, 1080) (SGN(I), I = 1, NCOOR)
         WRITE (LUPRI, 1090) (NR34H(I), I = 1, NCOOR)
         WRITE (LUPRI, 1100) IHC00
         WRITE (LUPRI, 1110) IHCHX, IHCHY, IHCHZ
         WRITE (LUPRI, 1120) NHC00
         WRITE (LUPRI, 1130) NHCHX, NHCHY, NHCHZ
      END IF
      RETURN
 1000 FORMAT(//,1X,'***** SUBROUTINE C1HINT - INITIALIZATION *****',/)
 1010 FORMAT(1X,'IEXADR  ',3I7)
 1020 FORMAT(1X,'NCOOR   ',I7)
 1030 FORMAT(1X,'NWORKH  ',I7)
 1040 FORMAT(1X,'NHCGT0  ',I7)
 1050 FORMAT(1X,'IADDHH  ',4I7)
 1060 FORMAT(1X,'JODDIF  ',4I7)
 1070 FORMAT(1X,'IOFFHC  ',4I7)
 1080 FORMAT(1X,'SGN     ',4F5.2)
 1090 FORMAT(1X,'NR34H   ',4I7)
 1100 FORMAT(1X,'IHC00   ',I7)
 1110 FORMAT(1X,'IHCHX   ',3I7)
 1120 FORMAT(1X,'NHC00   ',I7)
 1130 FORMAT(1X,'NHCHX   ',3I7)
      END
