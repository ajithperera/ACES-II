      SUBROUTINE GIACOE(BIGVEC,ONECEN,SEGMEN,
     &                  NHKTA,NHKTB,ISTEPA,ISTEPB,NUCAB,
     &                  IAB0X,IAB0Y,IAB0Z,DIFABX,DIFABY,DIFABZ,
     &                  SIGNAX,SIGNAY,SIGNAZ,SIGNBX,SIGNBY,SIGNBZ,
     &                  THRESH,ITYPE,MAXDER,DERECC,DERECD,DTEST,IPRINT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (ONE = 1.00 D00, D100 = 100.0 D00)
      LOGICAL YNEW, ZNEW, YEQX, ZEQX, ZEQY, ONECEN, DTEST,
     &        DERA, DERB, DERAB, DERECC, DERECD, BIGVEC, SEGMEN
#include "lwork.par"
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /TWOVEC/ DIFPAX(MXAOSQ), DIFPAY(MXAOSQ), DIFPAZ(MXAOSQ),
     &                DIFPBX(MXAOSQ), DIFPBY(MXAOSQ), DIFPBZ(MXAOSQ),
     &                TEXP1(MXAOSQ),  TEXP2(MXAOSQ),  HEXPPI(MXAOSQ)
      COMMON /ODCADR/ IE1X00, IE1X10, IE1X01, IE1X20, IE1X11, IE1X02,
     &                IE1Y00, IE1Y10, IE1Y01, IE1Y20, IE1Y11, IE1Y02,
     &                IE1Z00, IE1Z10, IE1Z01, IE1Z20, IE1Z11, IE1Z02,
     &                IE2X00, IE2X10, IE2X01, IE2X20, IE2X11, IE2X02,
     &                IE2Y00, IE2Y10, IE2Y01, IE2Y20, IE2Y11, IE2Y02,
     &                IE2Z00, IE2Z10, IE2Z01, IE2Z20, IE2Z11, IE2Z02,
     &                IELAST, ILST12
      COMMON /ODCCOM/ COORPX(MXAOSQ),COORQX(MXAOSQ),COORPY(MXAOSQ),
     &                COORQY(MXAOSQ),COORPZ(MXAOSQ),COORQZ(MXAOSQ),
     &                EXP12(MXAOSQ),EXP34(MXAOSQ),
     &                FAC12(MXAOSQ),FAC34(MXAOSQ)
      LA = NHKTA - 1
      LB = NHKTB - 1
      IF (ONECEN) THEN
       MAXDIF = 0
      ELSE
       MAXDIF = MAXDER
      END IF
      IF (ITYPE .EQ. 12) THEN
       SIGN = ONE
       IELAST = 0
       DERA = .TRUE.
       DERB = .TRUE.
       DERAB = .TRUE.
       IOFFC=1
      ELSE
       SIGN = - ONE
       IELAST = ILST12
       DERA = .TRUE.
       DERB = .TRUE.
       DERAB = .TRUE.
       IOFFC=1+MXAOSQ
      END IF
      MAXA = LA + MAXDIF
      MAXB = LB + MAXDIF
      MAXAB = LA + LB + MAXDIF
      ISTEPB = (MAXAB + 1)*NUCAB
      ISTEPA = (MAXB + 1)*ISTEPB
      ISTPAB = (MAXA + 1)*ISTEPA
C
C   DETERMINE UNIQUE COEFFICIENTS AND ADDRESSES 
C
C     ***** X-COEFFICIENTS *****
C
      IEX00 = IELAST
      IELAST = IELAST + ISTPAB
C
C     ***** Y-COEFFICIENTS *****
C
       IEY00 = IELAST
       IELAST = IELAST + ISTPAB
C
C     ***** Z-COEFFICIENTS *****
C
       IEZ00 = IELAST
       IELAST = IELAST + ISTPAB
C
C  CALCULATE UNDIFFERENTIATED EXPANSION COEFFICIENTS 
C
      IF (.NOT. (SEGMEN .AND. (MAXA+MAXB .EQ. 0)) ) THEN
       IADDAM = - ISTEPA - NUCAB
       IADDA0 = - ISTEPA
       IADDAP = - ISTEPA + NUCAB
       IADDBM = - ISTEPB - NUCAB
       IADDB0 = - ISTEPB
       IADDBP = - ISTEPB + NUCAB
       CALL TWOODC(MAXA,MAXB,IEX00,NUCAB,ISTEPA,ISTEPB,SIGN,SIGNAX,
     &             SIGNBX,IAB0X,IADDAM,IADDA0,IADDAP,IADDBM,IADDB0,
     &             IADDBP,DIFPAX,DIFPBX,HEXPPI,'EX00',IPRINT)
       CALL TWOODC(MAXA,MAXB,IEY00,NUCAB,ISTEPA,ISTEPB,SIGN,SIGNAY,
     &             SIGNBY,IAB0Y,IADDAM,IADDA0,IADDAP,IADDBM,IADDB0,
     &             IADDBP,DIFPAY,DIFPBY,HEXPPI,'EY00',IPRINT)
       CALL TWOODC(MAXA,MAXB,IEZ00,NUCAB,ISTEPA,ISTEPB,SIGN,SIGNAZ,
     &             SIGNBZ,IAB0Z,IADDAM,IADDA0,IADDAP,IADDBM,IADDB0,
     &             IADDBP,DIFPAZ,DIFPBZ,HEXPPI,'EZ00',IPRINT)
      END IF
C
C  CALCULATE FIRST DERIVATIVES OF OVERLAP DISTRIBUTIONS 
C
      IF (MAXDIF .GT. 0) THEN
       IEX10 = IELAST
       IELAST = IEX10 + ISTPAB
       CALL GIAONE(LA,LB,ISTEPA,ISTEPB,NUCAB,IEX10,IEX00,
     &             1,IAB0X,'EX10',IPRINT,COORPX(IOFFC),SIGN)
       IEY10 = IELAST
       IELAST = IEY10 + ISTPAB
       CALL GIAONE(LA,LB,ISTEPA,ISTEPB,NUCAB,IEY10,IEY00,
     &             1,IAB0Y,'EY10',IPRINT,COORPY(IOFFC),SIGN)
       IEZ10 = IELAST
       IELAST = IEZ10 + ISTPAB
       CALL GIAONE(LA,LB,ISTEPA,ISTEPB,NUCAB,IEZ10,IEZ00,
     &             1,IAB0Z,'EZ10',IPRINT,COORPZ(IOFFC),SIGN)
      END IF
C
C   CALCULATE SECOND DERIVATIVES OF OVERLAP DISTRIBUTIONS 
C
c      IF (MAXDIF .EQ. 2) THEN
c         IEX20 = IELAST
c         IEX02 = IEX20 + ISTPAB
c         IEX11 = IEX02 + ISTPAB
c         IELAST = IEX11 + ISTPAB
c         IF (DERA) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEX20,IEX00,
c     *                         1,IAB0X,'EX20',IPRINT)
c         IF (DERB) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEX02,IEX00,
c     *                         2,IAB0X,'EX02',IPRINT)
c         IF (DERAB) CALL DERTAB(LA,LB,ISTEPA,ISTEPB,NUCAB,IEX11,IEX00,
c     *                          IAB0X,'EX11',IPRINT)
c         IF (YEQX) THEN
c            IEY20 = IEX20
c            IEY02 = IEX02
c            IEY11 = IEX11
c         ELSE
c            IEY20 = IELAST
c            IEY02 = IEY20 + ISTPAB
c            IEY11 = IEY02 + ISTPAB
c            IELAST = IEY11 + ISTPAB
c            IF (DERA) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEY20,IEY00,
c     *                            1,IAB0Y,'EY20',IPRINT)
c            IF (DERB) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEY02,IEY00,
c     *                            2,IAB0Y,'EY02',IPRINT)
c            IF (DERAB) CALL DERTAB(LA,LB,ISTEPA,ISTEPB,NUCAB,IEY11,
c     *                             IEY00,IAB0Y,'EY11',IPRINT)
c         END IF
c         IF (ZEQX) THEN
c            IEZ20 = IEX20
c            IEZ02 = IEX02
c            IEZ11 = IEX11
c         ELSE IF (ZEQY) THEN
c            IEZ20 = IEY20
c            IEZ02 = IEY02
c            IEZ11 = IEY11
c         ELSE
c            IEZ20 = IELAST
c            IEZ02 = IEZ20 + ISTPAB
c            IEZ11 = IEZ02 + ISTPAB
c            IELAST = IEZ11 + ISTPAB
c            IF (DERA) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEZ20,IEZ00,
c     *                            1,IAB0Z,'EZ20',IPRINT)
c            IF (DERB) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEZ02,IEZ00,
c     *                            2,IAB0Z,'EZ02',IPRINT)
c            IF (DERAB) CALL DERTAB(LA,LB,ISTEPA,ISTEPB,NUCAB,IEZ11,
c     *                             IEZ00,IAB0Z,'EZ11',IPRINT)
c         END IF
c      END IF
C
C  SET WORK SPACE ADDRESSES 
C
      IF (ITYPE .EQ. 12) THEN
         IE1X00 = IEX00
         IE1Y00 = IEY00
         IE1Z00 = IEZ00
         IF (MAXDIF .GT. 0) THEN
            IE1X10 = IEX10
            IE1Y10 = IEY10
            IE1Z10 = IEZ10
c            IE1X01 = IEX01
c            IE1Y01 = IEY01
c            IE1Z01 = IEZ01
         END IF
c         IF (MAXDIF .EQ. 2) THEN
c            IE1X20 = IEX20
c            IE1Y20 = IEY20
c            IE1Z20 = IEZ20
c            IE1X02 = IEX02
c            IE1Y02 = IEY02
c            IE1Z02 = IEZ02
c            IE1X11 = IEX11
c            IE1Y11 = IEY11
c            IE1Z11 = IEZ11
c         END IF
         IF (IELAST .GT. LWORK2) THEN
            WRITE (LUPRI, 5000) IELAST, LWORK2
            STOP 
         END IF
         ILST12 = IELAST
      ELSE
         IE2X00 = IEX00
         IE2Y00 = IEY00
         IE2Z00 = IEZ00
         IF (MAXDIF .GT. 0) THEN
            IE2X10 = IEX10
            IE2Y10 = IEY10
            IE2Z10 = IEZ10
c            IE2X01 = IEX01
c            IE2Y01 = IEY01
c            IE2Z01 = IEZ01
         END IF
c         IF (MAXDIF .EQ. 2) THEN
c            IE2X20 = IEX20
c            IE2Y20 = IEY20
c            IE2Z20 = IEZ20
c            IE2X02 = IEX02
c            IE2Y02 = IEY02
c            IE2Z02 = IEZ02
c            IE2X11 = IEX11
c            IE2Y11 = IEY11
c            IE2Z11 = IEZ11
c         END IF
         IF (IELAST .GT. LWORK2) THEN
            WRITE (LUPRI, 5000) IELAST, LWORK2
            STOP 
         END IF
      END IF
      RETURN
 1000 FORMAT(//,' <<<<<<<<<< SUBROUTINE GIACOE >>>>>>>>>>',
     *       //,' Overlap distribution coefficients for electron ',I1,/)
 1010 FORMAT(/,'  NHKT: ',2I7)
 1020 FORMAT(/,'  SIGNA:',3I7, /, '  SIGNB:',3I7,/)
 1030 FORMAT('  ISTEPA ',I3,' ISTEPB ',I3,' ISTPAB ',I3,' NUCAB  ',I3)
 1040 FORMAT('  ITYPE: ',I3,' SIGN   ',F4.1)
 1050 FORMAT('  YNEW:  ',L3)
 1060 FORMAT('  ZNEW:  ',L3)
 1070 FORMAT('  YEQX:  ',L3)
 1080 FORMAT('  ZEQX:  ',L3)
 1090 FORMAT('  ZEQY:  ',L3)
 5000 FORMAT (//,'  work space requirement ',I6,' exceeds ',
     *        ' current limit ',I6,' of WORK2.')
      END
