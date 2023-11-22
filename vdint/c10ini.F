      SUBROUTINE C10INI
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
      LOGICAL HCGT0(1000), PATH1, SEGCON, TPRI12, TCON12
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
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     ***********************************************
C     ***** INITIALIZATION - ENTRY POINT C10INI *****
C     ***********************************************
C
      IF (PATH1) THEN
         IEXADR = IE1X00
         IEYADR = IE1Y00
         IEZADR = IE1Z00
      ELSE
         IEXADR = IE2X00
         IEYADR = IE2Y00
         IEZADR = IE2Z00
      END IF
      MAX340  = MAX34 + MAXDER
      NR34    = (MAX340 + 1)*(MAX340 + 2)*(MAX340 + 3)/6
      NRSTEP  = NR34*ISTEPT
      NRCCPP  = NR34*NCCPP
      IOFFHC  = IWKLST
      IDHC(1) = IWKLST
      NDHC(1) = NRCCPP
      IWKLST  = IWKLST + KHKT12*NRCCPP
      NWORK0  = NRSTEP
c      IF (IPRINT .GE. 10) THEN
c         WRITE (LUPRI, 1000)
c         WRITE (LUPRI, 1010) IEXADR, IEYADR, IEZADR
c         WRITE (LUPRI, 1020) NR34
c         WRITE (LUPRI, 1030) NRSTEP
c         WRITE (LUPRI, 1040) IOFFHC
c         WRITE (LUPRI, 1050) IWKLST
c      END IF
      RETURN
 1000 FORMAT(//,1X,'<<<<<<<<<< SUBROUTINE C10INT - INITIALIZATION ',
     *         '>>>>>>>>>>',/)
 1010 FORMAT(1X,'IEXADR ',3I7)
 1020 FORMAT(1X,'NR34   ',I7)
 1030 FORMAT(1X,'NRSTEP ',I7)
 1040 FORMAT(1X,'IOFFHC ',I7)
 1050 FORMAT(1X,'IWKLST ',I7)
      END
