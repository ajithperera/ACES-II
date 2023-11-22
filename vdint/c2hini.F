      SUBROUTINE C2HINI(WORK1,LWORK1)
C
C     TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (AZERO = 0.00 D00, ONE = 1.00 D00)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      LOGICAL PATH1, CROSS, DORDER(2), DHOVEC, DHTVEC, DMXVEC,
     *        PQSYM, DTEST, SEGCON, TPRI34, TCON34
      DIMENSION IXDER(6), IYDER(6), IZDER(6), JODD1(3),
     *          ITRI1(0:18), ITRI2(0:18),
     *          JODDIF(27), ITADD(27), IUADD(27), IVADD(27),
     *          FACTOR(27), IOFFHC(27), IOFFCC(27), NRCCPP(27),
     *          IDHO(2,2), IDHT(3,3,2), IDMX(3,3,2,2)
      DIMENSION WORK1(LWORK1)
#include "lwork.par"
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /CWORK3/ WK3LOW, WORK3(LWORK3), WK3HGH
      LOGICAL TKTIME
      COMMON /CC2INF/ KHKT12, KHKT34, KH1234,
     *                NORB3, NORB4, NORB12, NORB34, NO1234,
     *                TPRI34, TCON34,
     *                IPQ0X, IPQ0Y, IPQ0Z, INCRMT, INCRMU, INCRMV,
     *                ISTRET, ISTREU, ISTREV, LVAL34, MVAL34, NVAL34,
     *                NCCPP, NUC3, NUC4, NUC34, NSET3, NSET4,
     *                IADDPV, ISTRCF, IOFFVC, IODALL(MXAQN**2), PQSYM,
     *                ICMP34, IPATH, PATH1, MAXDER, CROSS, DTEST,
     *                SEGCON, IPRINT
      COMMON /PRIVEC/ NPRM12(MXAOSQ), NPRM34(MXAOSQ),
     *                LSTP12(MXAOSQ), LSTP34(MXAOSQ),
     *                LSTO12(MXAOSQ), LSTO34(MXAOSQ),
     *                LASTP(MXAOSQ),  LASTO(MXAOSQ)
      COMMON /INTADR/ IWKAO, IWKSO, IWKHHS, IWK1HH, IWK1HC, IWKLST
      COMMON /CCFCOM/ CONT1 (MXCONT*MXAOVC), CONT2 (MXCONT*MXAOVC),
     *                CONT3 (MXCONT*MXAOVC), CONT4 (MXCONT*MXAOVC),
     *                CONTT1(MXCONT*MXAOVC), CONTT2(MXCONT*MXAOVC),
     *                CONTT3(MXCONT*MXAOVC), CONTT4(MXCONT*MXAOVC),
     *                NUC1X (MXAOVC),        NUC2X (MXAOVC),
     *                NUC3X (MXAOVC),        NUC4X (MXAOVC),
     *                NRC1X (MXAOVC),        NRC2X (MXAOVC),
     *                NRC3X (MXAOVC),        NRC4X (MXAOVC)
      COMMON /ODCADR/ IE1X00, IE1X10, IE1X01, IE1X20, IE1X11, IE1X02,
     *                IE1Y00, IE1Y10, IE1Y01, IE1Y20, IE1Y11, IE1Y02,
     *                IE1Z00, IE1Z10, IE1Z01, IE1Z20, IE1Z11, IE1Z02,
     *                IE2X00, IE2X10, IE2X01, IE2X20, IE2X11, IE2X02,
     *                IE2Y00, IE2Y10, IE2Y01, IE2Y20, IE2Y11, IE2Y02,
     *                IE2Z00, IE2Z10, IE2Z01, IE2Z20, IE2Z11, IE2Z02,
     *                IELAST, ILST12
      LOGICAL         DZER
      COMMON /DERZER/ FZERO, DZER, IZERO
      COMMON /DHCINF/ IDHC(10), NDHC(10)
      COMMON /DHODIR/ DHOVEC(18)
      COMMON /DHOADR/ IHOVEC(18)
      COMMON /DHOFAC/ FHOVEC(18)
      COMMON /DHTDIR/ DHTVEC(9)
      COMMON /DHTADR/ IHTVEC(9)
      COMMON /DHTFAC/ FHTVEC(9)
      COMMON /DMXDIR/ DMXVEC(36)
      COMMON /DMXADR/ IMXVEC(36)
      COMMON /DMXFAC/ FMXVEC(36)
      COMMON/INT2H/NCOOR,JODDIF,ITADD,IUADD,IVADD,IOFFCC,FACTOR,
     *     IEXADR, IEYADR, IEZADR, NWORK, IADDCC, NRCCPP,
     *     IOFFHC
      DATA ITRI1 /0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,
     *            120,136,153,171/
     *     ITRI2 /0,1,4,10,20,35,56,84,120,165,220,286,364,455,
     *            560,680,816,969,1140/
      DATA IXDER /0,1,1,2,2,2/
     *     IYDER /0,1,0,2,1,0/
     *     IZDER /0,0,1,0,1,2/
      IHCADR(I,J,K) = NCCPP*(ITRI2(I+J+K) + ITRI1(J+K) + K)
C
C                *********************
C                ***** HOVEC(18) *****
C                *********************
C
      DATA IDHO /3,0,  12,6/
C
C     ARRANGEMENT OF VECTORS HOVEC(18)
C
C     1   XP00 YP00 ZP00
C     4   XQ00 YQ00 ZQ00
C     7   XXPP XYPP XZPP YYPP YZPP ZZPP
C     13  XXQQ XYQQ XZQQ YYQQ YZQQ ZZQQ
C
C                ********************
C                ***** HTVEC(9) *****
C                ********************
C
      DATA IDHT / 1, 2, 3, 4, 5, 6, 7, 8, 9,
     *            1, 4, 7, 2, 5, 8, 3, 6, 9/
C
C     ARRANGEMENT OF VECTORS HTVEC(9)
C
C     1   XXPQ XYPQ XZPQ
C     4   XYQP YYPQ YZPQ
C     7   XZQP YZQP ZZPQ
C
C                *********************
C                ***** MXVEC(36) *****
C                *********************
C
      DATA IDMX /19, 22, 25, 20, 23, 26, 21, 24, 27,
     *           28, 31, 34, 29, 32, 35, 30, 33, 36,
     *            1,  2,  3,  4,  5,  6,  7,  8,  9,
     *           10, 11, 12, 13, 14, 15, 16, 17, 18/
C
C     ARRANGEMENT OF VECTORS MXVEC(36)
C
C     1   XXPC XYPC XZPC
C     4   XYCP YYPC YZPC
C     7   XZCP YZCP ZZPC
C
C     10  XXPD XYPD XZPD
C     13  XYDP YYPD YZPD
C     16  XZDP YZDP ZZPD
C
C     19  XXAQ XYAQ XZAQ
C     22  XYQA YYAQ YZAQ
C     25  XZQA YZQA ZZAQ
C
C     28  XXBQ XYBQ XZBQ
C     31  XYQB YYBQ YZBQ
C     34  XZQB YZQB ZZBQ
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     *****************************************
C     ***** INITIALIZATION - ENTRY C2HINI *****
C     *****************************************
C
C     IEXADR, IEYADR & IEZADR
C
      IF (PATH1) THEN
         IEXADR = IE2X00
         IEYADR = IE2Y00
         IEZADR = IE2Z00
      ELSE
         IEXADR = IE1X00
         IEYADR = IE1Y00
         IEZADR = IE1Z00
      END IF
C
C     JODD1(3)
C
      IF (CROSS) THEN
         JODD1(1) = IPQ0X
         JODD1(2) = IBTSHL(IPQ0Y,1)
         JODD1(3) = IBTSHL(IPQ0Z,2)
      END IF
C
C     NCOOR, ITADD, IUADD, IVADD, JODDIF, IOFFCC & FACTOR
C
      ICOOR = 0
      IF (DZER) THEN
         ICOOR = ICOOR + 1
         ITADD(ICOOR) = 0
         IUADD(ICOOR) = 0
         IVADD(ICOOR) = 0
         NRCCPP(ICOOR) = NDHC(1)
         IOFFHC(ICOOR) = IDHC(1)
         JODDIF(ICOOR) = 0
         IOFFCC(ICOOR) = IWKAO + IZERO
         FACTOR(ICOOR) = FZERO
      END IF
      DORDER(1) = MAXDER .GE. 1
      DORDER(2) = MAXDER .EQ. 2
      IF (MAXDER .GT. 0) THEN
         DO 100 IORDER = 1, MAXDER
            IF (DORDER(IORDER)) THEN
               IF (PATH1 .AND. (IORDER .EQ. 1)) THEN
                  SGN = - ONE
               ELSE
                  SGN = ONE
               END IF
               JDHO = IDHO(IPATH,IORDER)
               DO 110 ICOOR2 = 1, (IORDER + 1)*(IORDER + 2)/2
                  IX = IORDER - IXDER(ICOOR2)
                  IY = IYDER(ICOOR2)
                  IZ = IZDER(ICOOR2)
                  JPQ0X = IBTAND(IPQ0X,IX)
                  JPQ0Y = IBTSHL(IBTAND(IPQ0Y,IY),1)
                  JPQ0Z = IBTSHL(IBTAND(IPQ0Z,IZ),2)
                  JODD2 = IBTOR(IBTOR(JPQ0Z,JPQ0Y),JPQ0X)
                  IADRHO = JDHO + ICOOR2
C
C                 DIFFERENTIATION ON A SINGLE ELECTRON
C
                  IADRHC = 1
                  IF (DHOVEC(IADRHO)) THEN
                     ICOOR = ICOOR + 1
                     ITADD(ICOOR) = IX
                     IUADD(ICOOR) = IY
                     IVADD(ICOOR) = IZ
                     NRCCPP(ICOOR) = NDHC(IADRHC)
                     IOFFHC(ICOOR) = IDHC(IADRHC)
                     JODDIF(ICOOR) = JODD2
                     IOFFCC(ICOOR) = IWKAO + IHOVEC(IADRHO)
                     FACTOR(ICOOR) = SGN*FHOVEC(IADRHO)
                  END IF
C
C                 DIFFERENTIATION ON TWO ELECTRONS
C
                  IF (DTEST) THEN
                     IF (CROSS .AND. (IORDER .EQ. 1)) THEN
                        DO 120 ICOOR1 = 1, 3
                           JODD12 = IBTXOR(JODD1(ICOOR1),JODD2)
                           IADRHT = IDHT(ICOOR2,ICOOR1,IPATH)
                           IADRHC = IADRHC + 1
                           IF (DHTVEC(IADRHT)) THEN
                              ICOOR = ICOOR + 1
                              ITADD(ICOOR) = IX
                              IUADD(ICOOR) = IY
                              IVADD(ICOOR) = IZ
                              NRCCPP(ICOOR) = NDHC(IADRHC)
                              IOFFHC(ICOOR) = IDHC(IADRHC)
                              JODDIF(ICOOR) = JODD12
                              IOFFCC(ICOOR) = IWKAO + IHTVEC(IADRHT)
                              FACTOR(ICOOR) = SGN*FHTVEC(IADRHT)
                           END IF
                           DO 130 IATOM1 = 1, 2
                              IADRMX = IDMX(ICOOR1,ICOOR2,IATOM1,IPATH)
                              IADRHC = IADRHC + 1
                              IF (DMXVEC(IADRMX)) THEN
                                 ICOOR = ICOOR + 1
                                 ITADD(ICOOR) = IX
                                 IUADD(ICOOR) = IY
                                 IVADD(ICOOR) = IZ
                                 NRCCPP(ICOOR) = NDHC(IADRHC)
                                 IOFFHC(ICOOR) = IDHC(IADRHC)
                                 JODDIF(ICOOR) = JODD12
                                 IOFFCC(ICOOR) = IWKAO + IMXVEC(IADRMX)
                                 FACTOR(ICOOR) = SGN*FMXVEC(IADRMX)
                              END IF
  130                      CONTINUE
  120                   CONTINUE
                     END IF
                  END IF
  110          CONTINUE
            END IF
  100    CONTINUE
      END IF
      NCOOR = ICOOR
      NWORK = NCOOR*KHKT12*NCCPP
      IWKMAX = IWKLST + NWORK
      IF (IWKMAX .GT. LWORK1) THEN
       WRITE (LUPRI, 1000) IWKMAX, IWKLST
       CALL ERREX
      END IF
      IADDCC = NORB12*KH1234
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 2000)
         WRITE (LUPRI, 2010)
         WRITE (LUPRI, 2020) DHOVEC(1)
         WRITE (LUPRI, 2030) (DHOVEC(I), I = 1, 3)
         WRITE (LUPRI, 2031) (DHOVEC(I), I = 4, 6)
         WRITE (LUPRI, 2032) (DHOVEC(I), I = 7, 12)
         WRITE (LUPRI, 2033) (DHOVEC(I), I = 13, 18)
         WRITE (LUPRI, 2034) (DHTVEC(I), I = 1, 9)
         WRITE (LUPRI, 2050) PATH1
         WRITE (LUPRI, 2060) KHKT12, KHKT34, KH1234
         WRITE (LUPRI, 2070) NORB12, NORB34
         WRITE (LUPRI, 2080) NUC34
         WRITE (LUPRI, 2090) NCCPP
         WRITE (LUPRI, 2110) IADDPV
         WRITE (LUPRI, 2114) NWORK
         WRITE (LUPRI, 2115) IWKMAX
         WRITE (LUPRI, 2120) IPQ0X, IPQ0Y, IPQ0Z
         WRITE (LUPRI, 2130) DTEST
         WRITE (LUPRI, 2140) NCOOR
         WRITE (LUPRI, 2150) (JODDIF(I), I = 1, NCOOR)
         WRITE (LUPRI, 2160) (ITADD(I), I = 1, NCOOR)
         WRITE (LUPRI, 2170) (IUADD(I), I = 1, NCOOR)
         WRITE (LUPRI, 2180) (IVADD(I), I = 1, NCOOR)
         WRITE (LUPRI, 2190) (IOFFCC(I), I = 1, NCOOR)
         WRITE (LUPRI, 2200) (FACTOR(I), I = 1, NCOOR)
         WRITE (LUPRI, 2210) (NRCCPP(I), I = 1, NCOOR)
      END IF
      RETURN
 1000 FORMAT (//,1X,' Work space requirement ',I8,' exceeds ',
     *        ' current limit ',I8,' of WORK1.')
 2000 FORMAT (//,1X,'<<<<<<<<<< SUBROUTINE C2HINT >>>>>>>>>>',/)
 2010 FORMAT (1X,   '           INITIALIZATION ',/)
 2020 FORMAT(1X,'DH0000',L7)
 2030 FORMAT(1X,'DHXP00',3L7)
 2031 FORMAT(1X,'DHXQ00',3L7)
 2032 FORMAT(1X,'DHXXPP',6L7)
 2033 FORMAT(1X,'DHXXQQ',6L7)
 2034 FORMAT(1X,'DHXXPQ',9L7)
 2050 FORMAT(1X,'PATH1 ',L7)
 2060 FORMAT(1X,'KHKT  ',3I7)
 2070 FORMAT(1X,'NORB  ',2I7)
 2080 FORMAT(1X,'NUC34 ',I7)
 2090 FORMAT(1X,'NCCPP ',I7)
 2110 FORMAT(1X,'IADDPV',I7)
 2114 FORMAT(1X,'NWORK ',I7)
 2115 FORMAT(1X,'IWKMAX',I7)
 2120 FORMAT(1X,'IPQ0X ',3I7)
 2130 FORMAT(1X,'DTEST ',L7)
 2140 FORMAT(1X,'NCOOR ',I7)
 2150 FORMAT(1X,'JODDIF',(10I7))
 2160 FORMAT(1X,'ITADD ',(10I7))
 2170 FORMAT(1X,'IUADD ',(10I7))
 2180 FORMAT(1X,'IVADD ',(10I7))
 2190 FORMAT(1X,'IOFFCC',(10I7))
 2200 FORMAT(1X,'FACTOR',(10F5.2))
 2210 FORMAT(1X,'NRCCPP',(10I7))
      END
