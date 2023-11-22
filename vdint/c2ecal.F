      SUBROUTINE C2ECAL(WORK1,LWORK1, Ienter)
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
      DO 100 ICOOR = 1, NCOOR
         MAXT = LVAL34 + INCMXT(ICOOR)
         MAXU = MVAL34 + INCMXU(ICOOR)
         MAXV = NVAL34 + INCMXV(ICOOR)
         MINT = IBTAND(MAXT,INCRMT - 1)
         MINU = IBTAND(MAXU,INCRMU - 1)
         MINV = IBTAND(MAXV,INCRMV - 1)
         ISTRAT = JSTRAT(ICOOR)
         IENDAT = JENDAT(ICOOR)
         NHCCMP = 0
         DO 110 ICMP12 = 1, KHKT12
            JODALL = IODALL(ICMP12)
            DO 120 IHCINT = JSTRHC(ISTRAT,ICOOR),JENDHC(IENDAT,ICOOR)
               IF (JODALL .EQ. JODDIF(IHCINT)) NHCCMP = NHCCMP + 1
  120       CONTINUE
  110    CONTINUE
         IF (NHCCMP .EQ. 0) GO TO 100
         NWORKP = NHCCMP*NCCPP
         NWORKC = NHCCMP*NO1234
C
C        *****************************************
C        ***** Primitive Cartesian Integrals *****
C        *****************************************
C
         CALL ZERO(WORK1(IWKLST+1),NWORKP)
         DO 200 IV = MINV, MAXV, INCRMV
            JSTREV = ISTREV + IV*NUC34
            DO 210 IU = MINU, MAXU, INCRMU
               JSTREU = ISTREU + IU*NUC34
               IUV = IU + IV
               IADRU = IV + ITRI1(IUV)
               DO 220 IT = MINT, MAXT, INCRMT
                  JSTRET = ISTRET + IT*NUC34
                  INTHC0 = NCCPP*(IADRU + ITRI2(IT + IUV))
                  IWORK = IWKLST
                  DO 300 IATOM = ISTRAT, IENDAT
                     ISTRHC = JSTRHC(IATOM,ICOOR)
                     IENDHC = JENDHC(IATOM,ICOOR)
                     KSTRET = IOFFEX(IATOM) + JSTRET
                     KSTREU = IOFFEY(IATOM) + JSTREU
                     KSTREV = IOFFEZ(IATOM) + JSTREV
                     DO 310 I = 1, NUC34
                        ECOEF(I) = WORK2(KSTRET + I)
     *                           * WORK2(KSTREU + I)
     *                           * WORK2(KSTREV + I)
  310                CONTINUE
C
C                    ***** KHKT12 = 1 *****
C
                     IF ((KHKT12 .EQ. 1) .AND. (ISTRHC.EQ.IENDHC)) THEN
                        IF (IODALL(1) .EQ. JODDIF(ISTRHC)) THEN
                           INTHC = IOFFHC(ISTRHC) + INTHC0
                           DO 320 I = 1, NORB12
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                              DO 330 J = 1, NUC34
                                 WORK1(IWORK + J) = WORK1(IWORK + J)
     *                                   + ECOEF(J)*WORK1(INTHC + J)
CSS                                 if (ienter .eq. 2) 
CSS     *                           Write(6,"(F12.6)") WORK1(IWORK + J) 
  330                         CONTINUE
                              INTHC = INTHC + NUC34
                              IWORK = IWORK + NUC34
  320                      CONTINUE
                        END IF
C
C                    ***** KHKT12 > 1 *****
C
                     ELSE
                        IJ = 0
                        DO 340 I = 1, NORB12
                           DO 350 J = 1, NUC34
                              IJ = IJ + 1
                              WORK3(IJ) = ECOEF(J)
  350                      CONTINUE
  340                   CONTINUE
                        DO 400 IHCINT = ISTRHC, IENDHC
                           INTHC  = IOFFHC(IHCINT) + INTHC0
                           IADHC0 = NRCCPP(IHCINT)
                           IODDIF = JODDIF(IHCINT)
                           DO 410 I = 1,KHKT12
                              IF (IODALL(I) .EQ. IODDIF) THEN
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                                 DO 420 J = 1, NCCPP
                                    WORK1(IWORK+J) = WORK1(IWORK+J)
     *                                  + WORK3(J)*WORK1(INTHC+J)
CSSS                                 if (ienter .eq. 2) 
CSSS     *                           Write(6,"(F12.6)") WORK1(IWORK + J) 
  420                            CONTINUE
                                 IWORK = IWORK + NCCPP
                              END IF
                              INTHC = INTHC + IADHC0
  410                      CONTINUE
  400                   CONTINUE
                     END IF
  300             CONTINUE
  220          CONTINUE
  210       CONTINUE
  200    CONTINUE
C
C        ******************************************
C        ***** Contracted Cartesian Integrals *****
C        ******************************************
C
         IWORK0 = IWKLST + 1
         INTCC0 = IWKLST + NWORKP + 1
         CALL ZERO(WORK1(INTCC0),NWORKC)
         KEND = NHCCMP*NORB12

CSSS         If (ienter .eq. 2) then 
CSSS            Write(6,*) "Primtive batch"
CSSS            write(6,"(6(1x,E12.6))") (WORK1(IWORK0+i),i=1,NSET4)
CSSS         endif

         CALL CNTINN(WORK1(IWORK0),WORK1(INTCC0),
     *               ISTRCF,IOFFVC,WORK3(1),
     *               NSET4,NSET3,NUC4,NUC3,NORB4,NORB3,
     *               KEND,TPRI34,TCON34)
CSSS         If (ienter .eq. 2) then 
CSSS            Write(6,*) "contracted batch"
CSSS            write(6,"(6(1x,E12.6))") (WORK1(INTCC0+i),i=1,NSET3)
CSSS         endif
C
C        Multiply by appropriate factors and distribute
C
C        Switch indices 12 and 34 for path 2
C
         IWORK = IWKLST + NWORKP
CSS         Write(6,*) PATH1
         IF (PATH1) THEN
            DO 800 IATOM = ISTRAT, IENDAT
               ISTRHC = JSTRHC(IATOM,ICOOR)
               IENDHC = JENDHC(IATOM,ICOOR)
               DO 810 IHCINT = ISTRHC, IENDHC
                  INTCC  = IOFFCC(IHCINT) + (ICMP34 - 1)*NO1234
                  FAC    = FACTOR(IHCINT)
                  IODDIF = JODDIF(IHCINT)
                  DO 820 ICMP12 = 1, KHKT12
                     IF (IODALL(ICMP12) .EQ. IODDIF) THEN
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                        DO 830 J = 1, NO1234
                           WORK1(INTCC + J) = WORK1(INTCC + J)
     *                                  + FAC*WORK1(IWORK + J)
  830                   CONTINUE
                        IWORK = IWORK + NO1234
                     END IF
                     INTCC = INTCC + KHKT34*NO1234
  820             CONTINUE
  810          CONTINUE
  800       CONTINUE
         ELSE
            IADDOF = (ICMP34 - 1)*KHKT12*NO1234 + 1
            DO 805 IATOM = ISTRAT, IENDAT
               ISTRHC = JSTRHC(IATOM,ICOOR)
               IENDHC = JENDHC(IATOM,ICOOR)
               DO 815 IHCINT = ISTRHC, IENDHC
                  ISTRCC = IOFFCC(IHCINT) + IADDOF
                  FAC    = FACTOR(IHCINT)
                  IODDIF = JODDIF(IHCINT)
                  DO 825 ICMP12 = 1, KHKT12
                  IF (IODALL(ICMP12) .EQ. IODDIF) THEN
                     INTCC0 = ISTRCC + (ICMP12 - 1)*NO1234
                     DO 835 I = 1, NORB12
                        INTCC = INTCC0
*VOCL LOOP,NOVREC
CDIR$ IVDEP
                        DO 845 J = 1, NORB34
                           WORK1(INTCC)=WORK1(INTCC)+FAC*WORK1(IWORK+J)
                           INTCC = INTCC + NORB12
  845                   CONTINUE
                        IWORK  = IWORK  + NORB34
                        INTCC0 = INTCC0 + 1
  835                CONTINUE
                  END IF
  825             CONTINUE
  815          CONTINUE
  805       CONTINUE
         END IF
         IPRINT = 0
         IF (IPRINT .GE. 10 .and. Ienter .EQ. 2) THEN
            WRITE (LUPRI,'(1X,A,3I7)') ' MINT  ', MINT, MINU, MINV
            WRITE (LUPRI,'(1X,A,3I7)') ' MAXT  ', MAXT, MAXU, MAXV
            WRITE (LUPRI,'(1X,A,3I7)') ' INCRMT', INCRMT,INCRMU,INCRMV
            IF (IPRINT .GE. 15) THEN
               CALL HEADER('Final Cartesian Integrals - C2EINT',-1)
               WRITE (LUPRI,'(/1X,A,I3)') ' ICMP34 ', ICMP34
               DO 900 IATOM = ISTRAT, IENDAT
                  WRITE (LUPRI,'(/1X,A,I3)') '  IATOM ', IATOM
                  DO 910 IHCINT=JSTRHC(IATOM,ICOOR),JENDHC(IATOM,ICOOR)
                     WRITE (LUPRI,'(/1X,A,I3)') ' IHCINT ', IHCINT
                     IF (PATH1) THEN
                        IADR = IOFFCC(IHCINT) + (ICMP34-1)*NO1234
                        IADD = KHKT34*NO1234
                     ELSE
                        IADR = IOFFCC(IHCINT) + (ICMP34-1)*KHKT12*NO1234
                        IADD = NO1234
                     END IF
                     DO 920 ICMP12 = 1, KHKT12
                        WRITE (LUPRI,'(/1X,A,I3)') ' ICMP12 ', ICMP12
                        IF (IODALL(ICMP12) .EQ. IODDIF) THEN
                           WRITE (LUPRI,'(1P,6E12.4)')
     *                           (WORK1(IADR+K),K=1,NO1234)
                        ELSE
                           WRITE (LUPRI,'(A)')
     *                        ' Integrals zero by symmetry'
                        END IF
                        IADR = IADR + IADD
  920                CONTINUE
  910             CONTINUE
  900          CONTINUE
            END IF
         END IF
  100 CONTINUE
      IPRINT = 0
      RETURN
 3000 FORMAT (//,1X,'<<<<<<<<<< SUBROUTINE C2EINT >>>>>>>>>>'/)
 3010 FORMAT (   1X,'           ENTRY POINT C2ECAL',/)
      END
