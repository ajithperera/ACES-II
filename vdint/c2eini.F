      SUBROUTINE C2EINI(WORK1,LWORK1,Ienter)
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
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      LOGICAL PATH1, CROSS, DCOOR, DCOOR0, DCOORH, DCOORE, DATOM,
     *        DATOM0, DATOMH, DATOME, DCORDH(2), DCORDE(2),
     *        DEOVEC, DETVEC, DMXVEC, DTEST,
     *        DEO, DET, DMX, PQSYM, SEGCON, TPRI34, TCON34
      DIMENSION ECOEF(MXAOSQ), ITRI1(0:18), ITRI2(0:18), IEOFF(0:2),
     *          IXDER(6), IYDER(6), IZDER(6), JODD1(3),
     *          INCMXT(9), INCMXU(9), INCMXV(9), JSTRAT(9),
     *          JENDAT(9), NWORK(9), JSTRHC(27,9), JENDHC(27,9),
     *          IOFFEX(27), IOFFEY(27), IOFFEZ(27),
     *          IOFFHC(81), IOFFCC(81), JODDIF(81), FACTOR(81),
     *          NRCCPP(81), IEND(6),
     *          IDEO(6,2,2), IDET(3,3,2,2,2), IDMX(3,3,2,2)
      DIMENSION WORK1(LWORK1)
#include "lwork.par"
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /CWORK3/ WK3LOW, WORK3(LWORK3), WK3HGH
      LOGICAL TKTIME
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
      COMMON /ODCADR/ IEADR(36), IELAST, ILST12
      COMMON /DEODIR/ DEOVEC(54)
      COMMON /DEOFAC/ FEOVEC(54)
      COMMON /DEOADR/ IEOVEC(54)
      COMMON /DETDIR/ DETVEC(36)
      COMMON /DETFAC/ FETVEC(36)
      COMMON /DETADR/ IETVEC(36)
      COMMON /DMXDIR/ DMXVEC(36)
      COMMON /DMXFAC/ FMXVEC(36)
      COMMON /DMXADR/ IMXVEC(36)
      COMMON /DHCINF/ IHCVEC(10), NHCVEC(10)
      COMMON /CC2INF/ KHKT12, KHKT34, KH1234,
     *                NORB3, NORB4, NORB12, NORB34, NO1234,
     *                TPRI34, TCON34,
     *                IPQ0X, IPQ0Y, IPQ0Z, INCRMT, INCRMU, INCRMV,
     *                ISTRET, ISTREU, ISTREV, LVAL34, MVAL34, NVAL34,
     *                NCCPP, NUC3, NUC4, NUC34, NSET3, NSET4,
     *                IADDPV, ISTRCF, IOFFVC, IODALL(MXAQN**2), PQSYM,
     *                ICMP34, IPATH, PATH1, MAXDER, CROSS, DTEST,
     *                SEGCON, IPRINT
      COMMON/INT2E/INCMXT,INCMXU,INCMXV,JSTRAT,JENDAT,JSTRHC,JENDHC,
     *             IOFFEX,IOFFEY,IOFFEZ,JODDIF,IOFFHC,IOFFCC,NRCCPP,
     *             FACTOR,NCOOR,NWORK
      DATA ITRI1 /0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,
     *            120,136,153,171/
      DATA ITRI2 /0,1,4,10,20,35,56,84,120,165,220,286,364,455,
     *            560,680,816,969,1140/
      DATA IXDER /0,1,1,2,2,2/
     *     IYDER /0,1,0,2,1,0/
     *     IZDER /0,0,1,0,1,2/
      DATA IEND  /3,4,4,3,4,3/
      DATA IEOFF /0,1,3/
C
C            ************************
C            ***** DEO(6,4,2,2) *****
C            ************************
C
      DEO(ICOOR2,IATOM2,IDER,IPATH)
     *     = DEOVEC(IDEO(ICOOR2,IDER,IPATH) + IATOM2)
C
      DATA IDEO / 2,  6, 10,  0,  0,  0,
     *           15, 22, 30, 37, 44, 51,
     *            0,  4,  8,  0,  0,  0,
     *           12, 18, 26, 34, 40, 48/
C
C     Arrangement of vector DEOVEC(54)
C
C     1   XA00 XB00             3   XC00 XD00
C     5   YA00 YB00             7   YC00 YD00
C     9   ZA00 ZB00             11  ZC00 ZD00
C     13  XXAA XXAB XXBB        16  XXCC XXCD XXDD
C     19  XYAA XYAB XYBA XYBB   23  XYCC XYCD XYDC XYDD
C     27  XZAA XZAB XZBA XZBB   31  XZCC XZCD XZDC XZDD
C     35  YYAA YYAB YYBB        38  YYCC YYCD YYDD
C     41  YZAA YZAB YZBA YZBB   45  YZCC YZCD YZDC YZDD
C     49  ZZAA ZZAB ZZBB        52  ZZCC ZZCD ZZDD
C
C            **************************
C            ***** DET(3,3,2,2,2) *****
C            **************************
C
      DET(ICOOR1,ICOOR2,IATOM1,IATOM2,IPATH)
     *     = DETVEC(IDET(ICOOR2,ICOOR1,IATOM2,IATOM1,IPATH))
C
      DATA IDET /  1,  2,  3,  4,  5,  6,  7,  8,  9,
     *            10, 11, 12, 13, 14, 15, 16, 17, 18,
     *            19, 20, 21, 22, 23, 24, 25, 26, 27,
     *            28, 29, 30, 31, 32, 33, 34, 35, 36,
     *             1,  4,  7,  2,  5,  8,  3,  6,  9,
     *            19, 22, 25, 20, 23, 26, 21, 24, 27,
     *            10, 13, 16, 11, 14, 17, 12, 15, 18,
     *            28, 31, 34, 29, 32, 35, 30, 33, 36/
C
C     Arrangement of vector DETVEC(36)
C
C     1   XXAC XYAC XZAC
C     4   XYCA YYAC YZAC
C     7   XZCA YZCA ZZAC
C
C     10  XXAD XYAD XZAD
C     13  XYDA YYAD YZAD
C     16  XZDA YZDA ZZAD
C
C     19  XXBC XYBC XZBC
C     22  XYCB YYBC YZBC
C     25  XZCB YZCB ZZBC
C
C     28  XXBD XYBD XZBD
C     31  XYDB YYBD YZBD
C     34  XZDB YZDB ZZBD
C
C            ************************
C            ***** DMX(3,3,2,2) *****
C            ************************
C
      DMX(ICOOR1,ICOOR2,IATOM2,IPATH)
     *     = DMXVEC(IDMX(ICOOR2,ICOOR1,IATOM2,IPATH))
C
      DATA IDMX / 1,  2,  3,  4,  5,  6,  7,  8,  9,
     *           10, 11, 12, 13, 14, 15, 16, 17, 18,
     *           19, 22, 25, 20, 23, 26, 21, 24, 27,
     *           28, 31, 34, 29, 32, 35, 30, 33, 36/
     *
C
C     Arrangement of vector DMXVEC(36)
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
C
C            **********************
C            ***** ODCADR(36) *****
C            **********************
C
C     Arrangement of elements in COMMON /ODCADR/
C
C     IE1X00 IE1X10 IE1X01 IE1X20 IE1X11 IE1X02
C     IE1Y00 IE1Y10 IE1Y01 IE1Y20 IE1Y11 IE1Y02
C     IE1Z00 IE1Z10 IE1Z01 IE1Z20 IE1Z11 IE1Z02
C     IE2X00 IE2X10 IE2X01 IE2X20 IE2X11 IE2X02
C     IE2Y00 IE2Y10 IE2Y01 IE2Y20 IE2Y11 IE2Y02
C     IE2Z00 IE2Y10 IE2Z01 IE2Z20 IE2Z11 IE2Z02
C
C
C            *********************
C            ***** HCVEC(10) *****
C            *********************
C
C     Arrangement of elements in COMMON /DHCINF/
C
C     HC00
C     HCHX HCEX1 HCEX2
C     HCHY HCEY1 HCEY2
C     HCHZ HCEZ1 HCEZ2
C
      IBTAND(I,J) = AND(I,J)
      IBTOR(I,J)  = OR(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IBTXOR(I,J) = XOR(I,J)
C
C     **************************
C     ***** Initialization *****
C     **************************
C
      IF (PATH1) THEN
         ISTRTX = 19
         ISTRTY = 25
         ISTRTZ = 31
      ELSE
         ISTRTX = 1
         ISTRTY = 7
         ISTRTZ = 13
      END IF
      IF (CROSS) THEN
         JODD1(1) = IPQ0X
         JODD1(2) = IBTSHL(IPQ0Y,1)
         JODD1(3) = IBTSHL(IPQ0Z,2)
      END IF
      IADDCC = NORB12*KH1234
      KHCCPP = KHKT12*NCCPP
C
C     MAX(ICOOR) = 3 + 6 = 9
C     MAX(IATOM) = 3*2 + (3 + 4 + 4 + 3 + 4 + 3) = 27
C     MAX(IHCINT) = 3*2*(1 + 3 + 3*2) + (3 + 4 + 4 + 3 + 4 + 3) = 81
C
      ICOOR = 0
      IATOM = 0
      IHCINT = 0
C
C     ***** First Derivatives and Cross Derivatives *****
C
      DO 10 ICOOR2 = 1, 3
         JDEO = IDEO(ICOOR2,1,IPATH)
         DCOOR0 = DEOVEC(JDEO + 1) .OR. DEOVEC(JDEO + 2)
         DCOOR = DCOOR0
         IF (CROSS) THEN
            DCORDH(1) = DMX(1,ICOOR2,1,IPATH)
     *             .OR. DMX(2,ICOOR2,1,IPATH)
     *             .OR. DMX(3,ICOOR2,1,IPATH)
            DCORDH(2) = DMX(1,ICOOR2,2,IPATH)
     *             .OR. DMX(2,ICOOR2,2,IPATH)
     *             .OR. DMX(3,ICOOR2,2,IPATH)
            DCORDE(1) = DET(1,ICOOR2,1,1,IPATH)
     *             .OR. DET(1,ICOOR2,2,1,IPATH)
     *             .OR. DET(2,ICOOR2,1,1,IPATH)
     *             .OR. DET(2,ICOOR2,2,1,IPATH)
     *             .OR. DET(3,ICOOR2,1,1,IPATH)
     *             .OR. DET(3,ICOOR2,2,1,IPATH)
            DCORDE(2) = DET(1,ICOOR2,1,2,IPATH)
     *             .OR. DET(1,ICOOR2,2,2,IPATH)
     *             .OR. DET(2,ICOOR2,1,2,IPATH)
     *             .OR. DET(2,ICOOR2,2,2,IPATH)
     *             .OR. DET(3,ICOOR2,1,2,IPATH)
     *             .OR. DET(3,ICOOR2,2,2,IPATH)
            DCOORH = DCORDH(1) .OR. DCORDH(2)
            DCOORE = DCORDE(1) .OR. DCORDE(2)
            DCOOR = DCOOR .OR. DCOORH .OR. DCOORE
         END IF
         IF (DCOOR) THEN
            ICOOR = ICOOR + 1
            JSTRAT(ICOOR) = IATOM + 1
            IWORK = IHCINT
            IX = 1 - IXDER(ICOOR2)
            IY = IYDER(ICOOR2)
            IZ = IZDER(ICOOR2)
            INCMXT(ICOOR) = IX
            INCMXU(ICOOR) = IY
            INCMXV(ICOOR) = IZ
            IEXOFF = ISTRTX + IEOFF(IX)
            IEYOFF = ISTRTY + IEOFF(IY)
            IEZOFF = ISTRTZ + IEOFF(IZ)
            JPQ0X = IBTAND(IPQ0X,IX)
            JPQ0Y = IBTSHL(IBTAND(IPQ0Y,IY),1)
            JPQ0Z = IBTSHL(IBTAND(IPQ0Z,IZ),2)
            JODD2 = IBTOR(IBTOR(JPQ0Z,JPQ0Y),JPQ0X)
            IATOM2 = 0
            DO 20 IATOMX = 0, IX
            DO 20 IATOMY = 0, IY
            DO 20 IATOMZ = 0, IZ
               IATOM2 = IATOM2 + 1
               IADREO = JDEO + IATOM2
               DATOM0 = DEOVEC(IADREO)
               DATOM = DATOM0
               IF (CROSS) THEN
                  DATOMH = DCORDH(IATOM2)
                  DATOME = DCORDE(IATOM2)
                  DATOM = DATOM .OR. DATOMH .OR. DATOME
               END IF
               IF (DATOM) THEN
                  IATOM = IATOM + 1
                  JSTRHC(IATOM,ICOOR) = IHCINT + 1
                  IOFFEX(IATOM) = IEADR(IEXOFF + IATOMX)
                  IOFFEY(IATOM) = IEADR(IEYOFF + IATOMY)
                  IOFFEZ(IATOM) = IEADR(IEZOFF + IATOMZ)
                  IF (DATOM0) THEN
                     IHCINT = IHCINT + 1
                     JODDIF(IHCINT) = JODD2
                     FACTOR(IHCINT) = FEOVEC(IADREO)
                     IOFFCC(IHCINT) = IWKAO + IEOVEC(IADREO)
                     IOFFHC(IHCINT) = IHCVEC(1)
                     NRCCPP(IHCINT) = NHCVEC(1)
                  END IF
                  IF (CROSS) THEN
                     DO 30 ICOOR1 = 1, 3
                        JODD = IBTXOR(JODD1(ICOOR1),JODD2)
                        IADRHC = 3*ICOOR1 - 1
                        IF (DATOMH) THEN
                           IADRMX = IDMX(ICOOR2,ICOOR1,IATOM2,IPATH)
                           IF (DMXVEC(IADRMX)) THEN
                              IHCINT = IHCINT + 1
                              JODDIF(IHCINT) = JODD
                              FACTOR(IHCINT) = FMXVEC(IADRMX)
                              IOFFCC(IHCINT) = IWKAO + IMXVEC(IADRMX)
                              IOFFHC(IHCINT) = IHCVEC(IADRHC)
                              NRCCPP(IHCINT) = NHCVEC(IADRHC)
                           END IF
                        END IF
                        IF (DATOME) THEN
                           DO 40 IATOM1 = 1, 2
                              IADRET =
     *                           IDET(ICOOR2,ICOOR1,IATOM2,IATOM1,IPATH)
                              IF (DETVEC(IADRET)) THEN
                                 IHCINT = IHCINT + 1
                                 JODDIF(IHCINT) = JODD
                                 FACTOR(IHCINT) = FETVEC(IADRET)
                                 IOFFCC(IHCINT) = IWKAO + IETVEC(IADRET)
                                 IOFFHC(IHCINT) = IHCVEC(IADRHC+IATOM1)
                                 NRCCPP(IHCINT) = NHCVEC(IADRHC+IATOM1)
                              END IF
   40                      CONTINUE
                        END IF
   30                CONTINUE
                  END IF
                 JENDHC(IATOM,ICOOR) = IHCINT
               END IF
   20       CONTINUE
            JENDAT(ICOOR) = IATOM
            NWORK(ICOOR) = KHCCPP*(IHCINT - IWORK)
         END IF
   10 CONTINUE
C
C     ***** Second Derivatives *****
C
      IF (MAXDER .EQ. 2) THEN
         DO 50 ICOOR2 = 1, 6
            JDEO = IDEO(ICOOR2,2,IPATH)
            DCOOR = .FALSE.
            DO 55 IATOM2 = 1, IEND(ICOOR2)
               DCOOR = DCOOR .OR. DEOVEC(JDEO + IATOM2)
   55       CONTINUE
            IF (DCOOR) THEN
               ICOOR = ICOOR + 1
               JSTRAT(ICOOR) = IATOM + 1
               IWORK = IHCINT
               IX = 2 - IXDER(ICOOR2)
               IY = IYDER(ICOOR2)
               IZ = IZDER(ICOOR2)
               INCMXT(ICOOR) = IX
               INCMXU(ICOOR) = IY
               INCMXV(ICOOR) = IZ
               IEXOFF = ISTRTX + IEOFF(IX)
               IEYOFF = ISTRTY + IEOFF(IY)
               IEZOFF = ISTRTZ + IEOFF(IZ)
               JPQ0X = IBTAND(IPQ0X,IX)
               JPQ0Y = IBTSHL(IBTAND(IPQ0Y,IY),1)
               JPQ0Z = IBTSHL(IBTAND(IPQ0Z,IZ),2)
               JODD2 = IBTOR(IBTOR(JPQ0Z,JPQ0Y),JPQ0X)
               IATOM2 = 0
               DO 60 IATOMX = 0, IX
               DO 60 IATOMY = 0, IY
               DO 60 IATOMZ = 0, IZ
                  IATOM2 = IATOM2 + 1
                  IADREO = JDEO + IATOM2
                  IF (DEOVEC(IADREO)) THEN
                     IATOM = IATOM + 1
                     IHCINT = IHCINT + 1
                     JSTRHC(IATOM,ICOOR) = IHCINT
                     JENDHC(IATOM,ICOOR) = IHCINT
                     IOFFEX(IATOM) = IEADR(IEXOFF + IATOMX)
                     IOFFEY(IATOM) = IEADR(IEYOFF + IATOMY)
                     IOFFEZ(IATOM) = IEADR(IEZOFF + IATOMZ)
                     JODDIF(IHCINT) = JODD2
                     FACTOR(IHCINT) = FEOVEC(IADREO)
                     IOFFCC(IHCINT) = IWKAO + IEOVEC(IADREO)
                     IOFFHC(IHCINT) = IHCVEC(1)
                     NRCCPP(IHCINT) = NHCVEC(1)
                  END IF
   60          CONTINUE
               JENDAT(ICOOR) = IATOM
               NWORK(ICOOR) = KHCCPP*(IHCINT - IWORK)
            END IF
   50    CONTINUE
      END IF
      NCOOR = ICOOR
      MAXWRK = 0
      DO 70 ICOOR = 1, NCOOR
         MAXWRK = MAX(MAXWRK,2*NWORK(ICOOR))
   70 CONTINUE
      IWKMAX = IWKLST + MAXWRK
      IF (IWKMAX .GT. LWORK1) THEN
       WRITE (LUPRI, 1000) IWKMAX, IWKLST
       CALL ERREX  
      END IF
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IPRINT = 0
      IF (IPRINT .GE. 10 .AND. Ienter .EQ. 2) THEN
         WRITE (LUPRI, 2000)
         WRITE (LUPRI, 2010)
         WRITE (LUPRI, 2020) PATH1
         WRITE (LUPRI, 2030) CROSS
         WRITE (LUPRI, 2040) MAXDER
         WRITE (LUPRI, 2050) KHKT12, KHKT34, KH1234
         WRITE (LUPRI, 2060) NORB12, NORB34
         WRITE (LUPRI, 2070) NCOOR
         WRITE (LUPRI, 2080) IATOM
         WRITE (LUPRI, 2090) IHCINT
         WRITE (LUPRI, 2095) MAXWRK
         WRITE (LUPRI, 2096) IWKMAX
         WRITE (LUPRI, 2100) (INCMXT(I), I = 1, NCOOR)
         WRITE (LUPRI, 2110) (INCMXU(I), I = 1, NCOOR)
         WRITE (LUPRI, 2120) (INCMXV(I), I = 1, NCOOR)
         WRITE (LUPRI, 2130) (JSTRAT(I), I = 1, NCOOR)
         WRITE (LUPRI, 2140) (JENDAT(I), I = 1, NCOOR)
         WRITE (LUPRI, 2145) (NWORK(I), I = 1, NCOOR)
         WRITE (LUPRI, 2150)
         DO 80 I = 1, NCOOR
            WRITE (LUPRI, 2151) (JSTRHC(J,I), J = JSTRAT(I),JENDAT(I))
   80    CONTINUE
         WRITE (LUPRI, 2160)
         DO 90 I = 1, NCOOR
            WRITE (LUPRI, 2161) (JENDHC(J,I), J = JSTRAT(I),JENDAT(I))
   90    CONTINUE
         WRITE (LUPRI, 2170) (IOFFEX(I), I = 1, IATOM)
         WRITE (LUPRI, 2180) (IOFFEY(I), I = 1, IATOM)
         WRITE (LUPRI, 2190) (IOFFEZ(I), I = 1, IATOM)
         WRITE (LUPRI, 2200) (IOFFHC(I), I = 1, IHCINT)
         WRITE (LUPRI, 2210) (IOFFCC(I), I = 1, IHCINT)
         WRITE (LUPRI, 2220) (JODDIF(I), I = 1, IHCINT)
         WRITE (LUPRI, 2230) (FACTOR(I), I = 1, IHCINT)
         WRITE (LUPRI, 2240) (NRCCPP(I), I = 1, IHCINT)
         WRITE (LUPRI, 2250) IADDCC
         WRITE (LUPRI, 2260) IPQ0X, IPQ0Y, IPQ0Z
         WRITE (LUPRI, 2270) INCRMT, INCRMU, INCRMV
         WRITE (LUPRI, 2280) NCCPP
         WRITE (LUPRI, 2290) NUC34
         WRITE (LUPRI, 2300) IADDPV
      END IF
      IPRINT = 0
      RETURN
 1000 FORMAT (//,1X,' Work space requirement ',I8,' exceeds ',
     *                      ' current limit ',I8,' of WORK1.')
 2000 FORMAT (//,1X,'<<<<<<<<<< SUBROUTINE C2EINT >>>>>>>>>>'/)
 2010 FORMAT (   1X,'           INITIALIZATION ',/)
 2020 FORMAT(1X,'PATH1 ',L7)
 2030 FORMAT(1X,'CROSS ',L7)
 2040 FORMAT(1X,'MAXDER',I7)
 2050 FORMAT(1X,'KHKT  ',3I7)
 2060 FORMAT(1X,'NORB  ',2I7)
 2070 FORMAT(1X,'NCOOR ',I7)
 2080 FORMAT(1X,'IATOM ',I7)
 2090 FORMAT(1X,'IHCINT',I7)
 2095 FORMAT(1X,'MAXWRK',I7)
 2096 FORMAT(1X,'IWKMAX',I7)
 2100 FORMAT(1X,'INCMXT',9(I7))
 2110 FORMAT(1X,'INCMXU',9(I7))
 2120 FORMAT(1X,'INCMXV',9(I7))
 2130 FORMAT(1X,'JSTRAT',9(I7))
 2140 FORMAT(1X,'JENDAT',9(I7))
 2145 FORMAT(1X,'NWORK ',9(I7))
 2150 FORMAT(/,1X,'JSTRHC')
 2151 FORMAT(7X,10I7)
 2160 FORMAT(/,1X,'JENDHC')
 2161 FORMAT(7X,10I7)
 2170 FORMAT(/,1X,'IOFFEX',(10I6))
 2180 FORMAT(1X,'IOFFEY',(10I6))
 2190 FORMAT(1X,'IOFFEZ',(10I6))
 2200 FORMAT(1X,'IOFFHC',(10I6))
 2210 FORMAT(1X,'IOFFCC',(10I6))
 2220 FORMAT(1X,'JODDIF',(10I6))
 2230 FORMAT(1X,'FACTOR',(10F6.2))
 2240 FORMAT(1X,'NRCCPP',(10I6))
 2250 FORMAT(1X,'IADDCC',I7)
 2260 FORMAT(1X,'IPQ0X ',3I7)
 2270 FORMAT(1X,'INCRMT',3I7)
 2280 FORMAT(1X,'NCCPP ',I7)
 2290 FORMAT(1X,'NUC34 ',I7)
 2300 FORMAT(1X,'IADDPV',I7)
      END
