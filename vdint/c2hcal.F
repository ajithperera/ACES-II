      SUBROUTINE C2HCAL(WORK1,LWORK1)
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
      NHCCMP = 0
      DO 200 ICMP12 = 1, KHKT12
         JODALL = IODALL(ICMP12)
         DO 210 ICOOR = 1, NCOOR
            IF (JODALL .EQ. JODDIF(ICOOR)) NHCCMP = NHCCMP + 1
  210    CONTINUE
  200 CONTINUE
      IF (NHCCMP .EQ. 0) RETURN
      NWORKP = NHCCMP*NCCPP
C
C     *******************************
C     ***** Primitive Integrals *****
C     *******************************
C
      CALL ZERO(WORK1(IWKLST+1),NWORKP)
      IF (LVAL34 + MVAL34 + NVAL34 .EQ. 0) THEN
         IWORK = IWKLST
         DO 605 ICOOR = 1, NCOOR
            IODDIF= JODDIF(ICOOR)
            JT = ITADD(ICOOR)
            JU = IUADD(ICOOR)
            JV = IVADD(ICOOR)
            INTHC0 = IOFFHC(ICOOR) + IHCADR(JT,JU,JV)
            IADHC0 = NRCCPP(ICOOR)
            DO 615 I = 1, KHKT12
               IF (IODALL(I) .EQ. IODDIF) THEN
                  INTHC = INTHC0 + (I - 1)*IADHC0
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                  DO 625 J = 1, NCCPP
                     WORK1(IWORK + J) = WORK1(INTHC + J)
  625             CONTINUE
                  IWORK = IWORK + NCCPP
               END IF
  615       CONTINUE
  605    CONTINUE
      ELSE
         MAXT = LVAL34
         MAXU = MVAL34
         MAXV = NVAL34
         MINT = IBTAND(MAXT,INCRMT - 1)
         MINU = IBTAND(MAXU,INCRMU - 1)
         MINV = IBTAND(MAXV,INCRMV - 1)
         IOFFEX = IEXADR + ISTRET
         IOFFEY = IEYADR + ISTREU
         IOFFEZ = IEZADR + ISTREV
         DO 400 IV = MINV, MAXV, INCRMV
            JSTREV = IOFFEZ + IV*NUC34
            DO 410 IU = MINU, MAXU, INCRMU
               JSTREU = IOFFEY + IU*NUC34
               DO 420 IT = MINT, MAXT, INCRMT
                  JSTRET = IOFFEX + IT*NUC34
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
                  DO 500 I = 1, NUC34
                     ECOEFI = WORK2(JSTRET + I)
     *                      * WORK2(JSTREU + I)
     *                      * WORK2(JSTREV + I)
                     IJ = I
                     DO 510 J = 1, NORB12
                        WORK3(IJ) = ECOEFI
                        IJ = IJ + NUC34
  510                CONTINUE
  500             CONTINUE
                  IWORK = IWKLST
                  DO 600 ICOOR = 1, NCOOR
                     IODDIF= JODDIF(ICOOR)
                     JT = IT + ITADD(ICOOR)
                     JU = IU + IUADD(ICOOR)
                     JV = IV + IVADD(ICOOR)
                     INTHC0 = IOFFHC(ICOOR) + IHCADR(JT,JU,JV)
                     IADHC0 = NRCCPP(ICOOR)
Cparallelization note: no synchronization problems
Cparallelization note: no synchronization problems
                     DO 610 I = 1, KHKT12
                        IF (IODALL(I) .EQ. IODDIF) THEN
                           INTHC = INTHC0 + (I - 1)*IADHC0
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                           DO 620 J = 1, NCCPP
                              WORK1(IWORK + J) = WORK1(IWORK + J)
     *                                + WORK3(J)*WORK1(INTHC + J)
  620                      CONTINUE
                           IWORK = IWORK + NCCPP
                        END IF
  610                CONTINUE
  600             CONTINUE
  420          CONTINUE
  410       CONTINUE
  400    CONTINUE
      END IF
C
C     ********************************
C     ***** Contracted Integrals *****
C     ********************************
C
      CALL CNTINN(WORK1(IWKLST+1),WORK1(IWKLST+NWORKP+1),
     *            ISTRCF,IOFFVC,WORK3(1),
     *            NSET4,NSET3,NUC4,NUC3,NORB4,NORB3,
     *            NHCCMP*NORB12,TPRI34,TCON34)
C
C     ***********************************************
C     ***** Multiply by factors and distribute ******
C     ***********************************************
C
      IWORK = IWKLST + NWORKP
      IF (PATH1) THEN
         DO 800 ICOOR = 1, NCOOR
            INTCC  = IOFFCC(ICOOR) + (ICMP34 - 1)*NO1234
            FAC    = FACTOR(ICOOR)
            IODDIF = JODDIF(ICOOR)
            DO 810 ICMP12 = 1, KHKT12
               IF (IODALL(ICMP12) .EQ. IODDIF) THEN
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                  DO 820 J = 1, NO1234
                     WORK1(INTCC + J) = WORK1(INTCC + J)
     *                            + FAC*WORK1(IWORK + J)
  820             CONTINUE
                  IWORK = IWORK + NO1234
               END IF
               INTCC = INTCC + KHKT34*NO1234
  810       CONTINUE
  800    CONTINUE
      ELSE
         IADDOF = (ICMP34 - 1)*KHKT12*NO1234 + 1
         DO 805 ICOOR = 1, NCOOR
            ISTRCC = IOFFCC(ICOOR) + IADDOF
            FAC    = FACTOR(ICOOR)
            IODDIF = JODDIF(ICOOR)
            DO 815 ICMP12 = 1, KHKT12
               IF (IODALL(ICMP12) .EQ. IODDIF) THEN
                  INTCC0 = ISTRCC + (ICMP12 - 1)*NO1234
                  DO 825 I = 1, NORB12
                     INTCC = INTCC0
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                     DO 835 J = 1, NORB34
                        WORK1(INTCC) = WORK1(INTCC) + FAC*WORK1(IWORK+J)
                        INTCC = INTCC + NORB12
  835                CONTINUE
                     INTCC0 = INTCC0 + 1
                     IWORK  = IWORK  + NORB34
  825             CONTINUE
               END IF
  815       CONTINUE
  805    CONTINUE
      END IF
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 3000)
         WRITE (LUPRI,'(1X,A,   I7 )') 'ICMP34', ICMP34
         WRITE (LUPRI,'(1X,A,  3I7 )') 'LVAL34', LVAL34, MVAL34, NVAL34
         WRITE (LUPRI,'(1X,A,  3I7 )') 'ISTRET', ISTRET, ISTREU, ISTREV
         WRITE (LUPRI,'(1X,A,  3I7 )') 'IOFFEX', IOFFEX, IOFFEY, IOFFEZ
         WRITE (LUPRI,'(1X,A,  3I7 )') 'MINT  ', MINT, MINU, MINV
         WRITE (LUPRI,'(1X,A,  3I7 )') 'MAXT  ', MAXT, MAXU, MAXV
         WRITE (LUPRI,'(1X,A,  3I7 )') 'INCRMT', INCRMT, INCRMU, INCRMV
         WRITE (LUPRI,'(1X,A,(10I7))') 'IODALL', (IODALL(I),I=1,KHKT12)
         IF (IPRINT .GE. 15) THEN
            CALL HEADER('Final Cartesian Integrals - C2HINT',-1)
            DO 900 ICOOR = 1, NCOOR
               WRITE (LUPRI, '(/1X,A,I3)') ' ICOOR ', ICOOR
               IF (PATH1) THEN
                  IADR = IOFFCC(ICOOR) + (ICMP34-1)*NO1234
                  IADD = KHKT34*NO1234
               ELSE
                  IADR = IOFFCC(ICOOR) + (ICMP34-1)*KHKT12*NO1234
                  IADD = NO1234
               END IF
               DO 910 ICMP12 = 1, KHKT12
                  WRITE (LUPRI,'(/1X,A,I3/)') ' ICMP12 ', ICMP12
                  IF (IODALL(ICMP12) .EQ. IODDIF) THEN
                    WRITE(LUPRI,'(1P,6E12.4)')(WORK1(IADR+K),K=1,NO1234)
                  ELSE
                     WRITE(LUPRI,'(A)') ' Integrals zero by symmetry'
                  END IF
                  IADR = IADR + IADD
  910          CONTINUE
  900       CONTINUE
         END IF
      END IF
      RETURN
 3000 FORMAT (//,1X,'<<<<<<<<<< SUBROUTINE C2HINT >>>>>>>>>>'/)
      END
