      SUBROUTINE C1ECAL(WORK1,LWORK1)
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
C     *********************************************************
C     ***** CALCULATION OF INTEGRALS - ENTRY POINT C1ECAL *****
C     *********************************************************
C
C
C     *****************************************
C     ***** PRIMITIVE CARTESIAN INTEGRALS *****
C     *****************************************
C
      DO 100 ICOOR = 1, NCOOR
         MAXT = LVAL12 + INCMXT(ICOOR)
         MAXU = MVAL12 + INCMXU(ICOOR)
         MAXV = NVAL12 + INCMXV(ICOOR)
         MINT = IBTAND(MAXT,INCRMT - 1)
         MINU = IBTAND(MAXU,INCRMU - 1)
         MINV = IBTAND(MAXV,INCRMV - 1)
         ISTRAT = JSTRAT(ICOOR)
         IENDAT = JENDAT(ICOOR)
c         IF (IPRINT .GE. 10) THEN
c            WRITE (LUPRI, 2010) ICOOR
c            WRITE (LUPRI, 2020) MINT, MAXT, INCRMT
c            WRITE (LUPRI, 2030) MINU, MAXU, INCRMU
c            WRITE (LUPRI, 2040) MINV, MAXV, INCRMV
c            WRITE (LUPRI, 2060) ISTRAT
c            WRITE (LUPRI, 2070) IENDAT
c            WRITE (LUPRI, 2080) NR34
c         END IF
C
C        *****************************************
C        ***** PRIMITIVE CARTESIAN INTEGRALS *****
C        *****************************************
C
         DO 150 I = 1, NR34
            HCGT0(I) = .FALSE.
  150    CONTINUE
         DO 200 IATOM = ISTRAT, IENDAT
            IOFFEX = JOFFEX(IATOM) + ISTRET
            IOFFEY = JOFFEY(IATOM) + ISTREU
            IOFFEZ = JOFFEZ(IATOM) + ISTREV
            DO 250 I = 1, NWORKE
               WORK1(IWKLST + I) = AZERO
  250       CONTINUE
            DO 300 IV = MINV, MAXV, INCRMV
               JSTREV = IOFFEZ + IV*NUC12
               ISTRHV = IV*ISTEPV
               IODDZ = IBTSHL(IBTAND(IPQ0Z,IV),2)
               DO 310 IU = MINU, MAXU, INCRMU
                  JSTREU = IOFFEY + IU*NUC12
                  ISTRHU = ISTRHV + IU*ISTEPU
                  IODDYZ = IBTOR(IODDZ,IBTSHL(IBTAND(IPQ0Y,IU),1))
                  DO 320 IT = MINT, MAXT, INCRMT
                     JSTRET = IOFFEX + IT*NUC12
                     INTHH0 = ISTRHU + IT
                     IODXYZ = IBTOR(IODDYZ,IBTAND(IPQ0X,IT))
                     IJ = 0
                     DO 400 I = 1, NUC12
                        ECOEFI = WORK2(JSTRET + I)
     *                         * WORK2(JSTREU + I)
     *                         * WORK2(JSTREV + I)
                        DO 410 J = 1, NUC34
                           IJ = IJ + 1
                           WORK3(IJ) = ECOEFI
  410                   CONTINUE
  400                CONTINUE
                     IWORK = IWKLST
                     DO 500 ITUV = 1, NR34
                        IF (JODDH(ITUV) .EQ. IODXYZ) THEN
                           HCGT0(ITUV) = .TRUE.
                           INTHH = INAHH(INTHH0 + JSTRH(ITUV))
CDIR$ IVDEP
*VOCL LOOP,NOVREC
                           DO 600 I = 1, ISTEPT
                              WORK1(IWORK + I) = WORK1(IWORK + I)
     *                                + WORK3(I)*WORK1(INTHH + I)
  600                      CONTINUE
                        END IF
                        IWORK = IWORK + ISTEPT
  500                CONTINUE
  320             CONTINUE
  310          CONTINUE
  300       CONTINUE
C
C           ******************************************
C           ***** CONTRACTED CARTESIAN INTEGRALS *****
C           ******************************************
C
            IWORK0 = IWKLST
            INTHC0 = IOFFHC(IATOM)
            CALL ZERO(WORK1(INTHC0 + 1),NRCCPP)
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
     *              ' Work space for CNTOUT insufficient.',
     *              ' Current value of LWORK3',LWORK3,
     *              ' must be increased to at least ',LSCR12
               STOP 
            END IF
            ISCR1 = 1
            ISCR2 = ISCR1 + LSCR2
            DO 700 ITUV = 1, NR34
               IF (HCGT0(ITUV)) THEN
                  CALL CNTOUT(WORK1(IWORK0+1),WORK1(INTHC0+1),
     *                        WORK3(ISCR1),WORK3(ISCR2),LSCR1,LSCR2,
     *                        ISTRCF,IOFFVC,
     *                        NSET2,NSET1,NUC2,NUC1,NORB2,NORB1,
     *                        NUC34,TPRI12,TCON12)
c                  IF (IPRINT .GE. 20) THEN
c                     WRITE (LUPRI, 3000) IATOM
c                     WRITE (LUPRI, 3010) ITUV
c                     WRITE (LUPRI, 3020) IWORK0
c                     WRITE (LUPRI, 3030) INTHC0
c                     INTHC = INTHC0
c                     DO 900 I = 1, NORB12
c                        WRITE (LUPRI, 3040) I,
c     *                        (WORK1(INTHC+J), J=1,NUC34)
c                        INTHC = INTHC + NUC34
c  900                CONTINUE
c                  END IF
               END IF
               INTHC0 = INTHC0 + NCCPP
               IWORK0 = IWORK0 + ISTEPT
  700       CONTINUE
            IOFFHC(IATOM) = INTHC0
  200    CONTINUE
  100 CONTINUE
      RETURN
 2000 FORMAT (/,1X,'***** EC-DIFFERENTIATED HC-INTEGRALS *****',/)
 2010 FORMAT (1X,'ICOOR  ',I7)
 2020 FORMAT (1X,'LOOP PARAMETERS X-DIR ',3I7)
 2030 FORMAT (1X,'LOOP PARAMETERS Y-DIR ',3I7)
 2040 FORMAT (1X,'LOOP PARAMETERS Z-DIR ',3I7)
 2060 FORMAT (1X,'ISTRAT ',I7)
 2070 FORMAT (1X,'IENDAT ',I7)
 2080 FORMAT (1X,'NR34   ',I7)
 3000 FORMAT (1X,'IATOM  ',I7)
 3010 FORMAT (1X,'ITUV   ',I7)
 3020 FORMAT (1X,'IWORK0 ',I7)
 3030 FORMAT (1X,'INTHC0 ',I7)
 3040 FORMAT (/,1X,'ORB12',I7/,(1P,6E12.4))
      END
