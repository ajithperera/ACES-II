      SUBROUTINE EXCOEF(BIGVEC,ONECEN,SEGMEN,
     *                  NHKTA,NHKTB,ISTEPA,ISTEPB,NUCAB,
     *                  IAB0X,IAB0Y,IAB0Z,DIFABX,DIFABY,DIFABZ,
     *                  SIGNAX,SIGNAY,SIGNAZ,SIGNBX,SIGNBY,SIGNBZ,
     *                  THRESH,ITYPE,MAXDER,DERECC,DERECD,DTEST,IPRINT)
C
C     TUH 84
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (ONE = 1.00 D00, D100 = 100.0 D00)
      LOGICAL YNEW, ZNEW, YEQX, ZEQX, ZEQY, ONECEN, DTEST,
     *        DERA, DERB, DERAB, DERECC, DERECD, BIGVEC, SEGMEN
#include "lwork.par"
      COMMON /CWORK2/ WK2LOW, WORK2(LWORK2), WK2HGH
      COMMON /TWOVEC/ DIFPAX(MXAOSQ), DIFPAY(MXAOSQ), DIFPAZ(MXAOSQ),
     *                DIFPBX(MXAOSQ), DIFPBY(MXAOSQ), DIFPBZ(MXAOSQ),
     *                TEXP1(MXAOSQ),  TEXP2(MXAOSQ),  HEXPPI(MXAOSQ)
      COMMON /ODCADR/ IE1X00, IE1X10, IE1X01, IE1X20, IE1X11, IE1X02,
     *                IE1Y00, IE1Y10, IE1Y01, IE1Y20, IE1Y11, IE1Y02,
     *                IE1Z00, IE1Z10, IE1Z01, IE1Z20, IE1Z11, IE1Z02,
     *                IE2X00, IE2X10, IE2X01, IE2X20, IE2X11, IE2X02,
     *                IE2Y00, IE2Y10, IE2Y01, IE2Y20, IE2Y11, IE2Y02,
     *                IE2Z00, IE2Z10, IE2Z01, IE2Z20, IE2Z11, IE2Z02,
     *                IELAST, ILST12
      LA = NHKTA - 1
      LB = NHKTB - 1
      IF (ONECEN .AND. (.NOT.DTEST)) THEN
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
      ELSE
         SIGN = - ONE
         IELAST = ILST12
         DERA = DERECC .OR. DTEST
         DERB = DERECD .OR. DTEST
         DERAB = DERA .AND. DERB
      END IF
      MAXA = LA + MAXDIF
      MAXB = LB + MAXDIF
      MAXAB = LA + LB + MAXDIF
      ISTEPB = (MAXAB + 1)*NUCAB
      ISTEPA = (MAXB + 1)*ISTEPB
      ISTPAB = (MAXA + 1)*ISTEPA
C
C     *******************************************************
C     ***** DETERMINE UNIQUE COEFFICIENTS AND ADDRESSES *****
C     *******************************************************
C
      IF (BIGVEC) THEN
         DIFXY = D100
         DIFXZ = D100
         DIFYZ = D100
      ELSE
         DIFXY = ABS(DIFABX - DIFABY)
         DIFXZ = ABS(DIFABX - DIFABZ)
         DIFYZ = ABS(DIFABY - DIFABZ)
      END IF
C
C     ***** X-COEFFICIENTS *****
C
      IEX00 = IELAST
      IELAST = IELAST + ISTPAB
C
C     ***** Y-COEFFICIENTS *****
C
      IF (DIFXY .LT. THRESH) THEN
         IEY00 = IEX00
         YEQX = .TRUE.
      ELSE
         IEY00 = IELAST
         IELAST = IELAST + ISTPAB
         YEQX = .FALSE.
      END IF
      YNEW = .NOT. YEQX
C
C     ***** Z-COEFFICIENTS *****
C
      IF (DIFXZ .LT. THRESH) THEN
         IEZ00 = IEX00
         ZEQX = .TRUE.
         ZEQY = .FALSE.
      ELSE IF (DIFYZ .LT. THRESH) THEN
         IEZ00 = IEY00
         ZEQX = .FALSE.
         ZEQY = .TRUE.
      ELSE
         IEZ00 = IELAST
         IELAST = IELAST + ISTPAB
         ZEQX = .FALSE.
         ZEQY = .FALSE.
      END IF
      ZNEW = .NOT. (ZEQX .OR. ZEQY)
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GE. 10) THEN
         IELECT = 1
         IF (ITYPE .EQ. 34) IELECT = 2
         WRITE (LUPRI, 1000) IELECT
         WRITE (LUPRI, 1010) NHKTA, NHKTB
         WRITE (LUPRI, 1020)IDINT(SIGNAX), IDINT(SIGNAY), IDINT(SIGNAZ),
     *                      IDINT(SIGNBX), IDINT(SIGNBY), IDINT(SIGNBZ)
         WRITE (LUPRI, 1030) ISTEPA, ISTEPB, ISTPAB, NUCAB
         WRITE (LUPRI, 1040) ITYPE, SIGN
         WRITE (LUPRI, 1050) YNEW
         WRITE (LUPRI, 1060) ZNEW
         WRITE (LUPRI, 1070) YEQX
         WRITE (LUPRI, 1080) ZEQX
         WRITE (LUPRI, 1090) ZEQY
      END IF
C
C     *************************************************************
C     ***** CALCULATE UNDIFFERENTIATED EXPANSION COEFFICIENTS *****
C     *************************************************************
C
      IF (.NOT. (SEGMEN .AND. (MAXA+MAXB .EQ. 0)) ) THEN
         IADDAM = - ISTEPA - NUCAB
         IADDA0 = - ISTEPA
         IADDAP = - ISTEPA + NUCAB
         IADDBM = - ISTEPB - NUCAB
         IADDB0 = - ISTEPB
         IADDBP = - ISTEPB + NUCAB
         CALL TWOODC(MAXA,MAXB,IEX00,NUCAB,ISTEPA,ISTEPB,SIGN,SIGNAX,
     *               SIGNBX,IAB0X,IADDAM,IADDA0,IADDAP,IADDBM,IADDB0,
     *               IADDBP,DIFPAX,DIFPBX,HEXPPI,'EX00',IPRINT)
         IF (YNEW) THEN
            CALL TWOODC(MAXA,MAXB,IEY00,NUCAB,ISTEPA,ISTEPB,SIGN,SIGNAY,
     *                  SIGNBY,IAB0Y,IADDAM,IADDA0,IADDAP,IADDBM,IADDB0,
     *                  IADDBP,DIFPAY,DIFPBY,HEXPPI,'EY00',IPRINT)
         END IF
         IF (ZNEW) THEN
            CALL TWOODC(MAXA,MAXB,IEZ00,NUCAB,ISTEPA,ISTEPB,SIGN,SIGNAZ,
     *                  SIGNBZ,IAB0Z,IADDAM,IADDA0,IADDAP,IADDBM,IADDB0,
     *                  IADDBP,DIFPAZ,DIFPBZ,HEXPPI,'EZ00',IPRINT)
         END IF
      END IF
C
C     ****************************************************************
C     ***** CALCULATE FIRST DERIVATIVES OF OVERLAP DISTRIBUTIONS *****
C     ****************************************************************
C
      IF (MAXDIF .GT. 0) THEN
         IEX10 = IELAST
         IEX01 = IEX10 + ISTPAB
         IELAST = IEX01 + ISTPAB
         IF (DERA) CALL DERONE(LA,LB,ISTEPA,ISTEPB,NUCAB,IEX10,IEX00,
     *                         1,IAB0X,'EX10',IPRINT)
         IF (DERB) CALL DERONE(LA,LB,ISTEPA,ISTEPB,NUCAB,IEX01,IEX00,
     *                         2,IAB0X,'EX01',IPRINT)
         IF (YEQX) THEN
            IEY10 = IEX10
            IEY01 = IEX01
         ELSE
            IEY10 = IELAST
            IEY01 = IEY10 + ISTPAB
            IELAST = IEY01 + ISTPAB
            IF (DERA) CALL DERONE(LA,LB,ISTEPA,ISTEPB,NUCAB,IEY10,IEY00,
     *                            1,IAB0Y,'EY10',IPRINT)
            IF (DERB) CALL DERONE(LA,LB,ISTEPA,ISTEPB,NUCAB,IEY01,IEY00,
     *                            2,IAB0Y,'EY01',IPRINT)
         END IF
         IF (ZEQX) THEN
            IEZ10 = IEX10
            IEZ01 = IEX01
         ELSE IF (ZEQY) THEN
            IEZ10 = IEY10
            IEZ01 = IEY01
         ELSE
            IEZ10 = IELAST
            IEZ01 = IEZ10 + ISTPAB
            IELAST = IEZ01 + ISTPAB
            IF (DERA) CALL DERONE(LA,LB,ISTEPA,ISTEPB,NUCAB,IEZ10,IEZ00,
     *                            1,IAB0Z,'EZ10',IPRINT)
            IF (DERB) CALL DERONE(LA,LB,ISTEPA,ISTEPB,NUCAB,IEZ01,IEZ00,
     *                            2,IAB0Z,'EZ01',IPRINT)
         END IF
      END IF
C
C     *****************************************************************
C     ***** CALCULATE SECOND DERIVATIVES OF OVERLAP DISTRIBUTIONS *****
C     *****************************************************************
C
      IF (MAXDIF .EQ. 2) THEN
         IEX20 = IELAST
         IEX02 = IEX20 + ISTPAB
         IEX11 = IEX02 + ISTPAB
         IELAST = IEX11 + ISTPAB
         IF (DERA) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEX20,IEX00,
     *                         1,IAB0X,'EX20',IPRINT)
         IF (DERB) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEX02,IEX00,
     *                         2,IAB0X,'EX02',IPRINT)
         IF (DERAB) CALL DERTAB(LA,LB,ISTEPA,ISTEPB,NUCAB,IEX11,IEX00,
     *                          IAB0X,'EX11',IPRINT)
         IF (YEQX) THEN
            IEY20 = IEX20
            IEY02 = IEX02
            IEY11 = IEX11
         ELSE
            IEY20 = IELAST
            IEY02 = IEY20 + ISTPAB
            IEY11 = IEY02 + ISTPAB
            IELAST = IEY11 + ISTPAB
            IF (DERA) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEY20,IEY00,
     *                            1,IAB0Y,'EY20',IPRINT)
            IF (DERB) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEY02,IEY00,
     *                            2,IAB0Y,'EY02',IPRINT)
            IF (DERAB) CALL DERTAB(LA,LB,ISTEPA,ISTEPB,NUCAB,IEY11,
     *                             IEY00,IAB0Y,'EY11',IPRINT)
         END IF
         IF (ZEQX) THEN
            IEZ20 = IEX20
            IEZ02 = IEX02
            IEZ11 = IEX11
         ELSE IF (ZEQY) THEN
            IEZ20 = IEY20
            IEZ02 = IEY02
            IEZ11 = IEY11
         ELSE
            IEZ20 = IELAST
            IEZ02 = IEZ20 + ISTPAB
            IEZ11 = IEZ02 + ISTPAB
            IELAST = IEZ11 + ISTPAB
            IF (DERA) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEZ20,IEZ00,
     *                            1,IAB0Z,'EZ20',IPRINT)
            IF (DERB) CALL DERTAA(LA,LB,ISTEPA,ISTEPB,NUCAB,IEZ02,IEZ00,
     *                            2,IAB0Z,'EZ02',IPRINT)
            IF (DERAB) CALL DERTAB(LA,LB,ISTEPA,ISTEPB,NUCAB,IEZ11,
     *                             IEZ00,IAB0Z,'EZ11',IPRINT)
         END IF
      END IF
C
C     ************************************
C     ***** SET WORK SPACE ADDRESSES *****
C     ************************************
C
      IF (ITYPE .EQ. 12) THEN
         IE1X00 = IEX00
         IE1Y00 = IEY00
         IE1Z00 = IEZ00
         IF (MAXDIF .GT. 0) THEN
            IE1X10 = IEX10
            IE1Y10 = IEY10
            IE1Z10 = IEZ10
            IE1X01 = IEX01
            IE1Y01 = IEY01
            IE1Z01 = IEZ01
         END IF
         IF (MAXDIF .EQ. 2) THEN
            IE1X20 = IEX20
            IE1Y20 = IEY20
            IE1Z20 = IEZ20
            IE1X02 = IEX02
            IE1Y02 = IEY02
            IE1Z02 = IEZ02
            IE1X11 = IEX11
            IE1Y11 = IEY11
            IE1Z11 = IEZ11
         END IF
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
            IE2X01 = IEX01
            IE2Y01 = IEY01
            IE2Z01 = IEZ01
         END IF
         IF (MAXDIF .EQ. 2) THEN
            IE2X20 = IEX20
            IE2Y20 = IEY20
            IE2Z20 = IEZ20
            IE2X02 = IEX02
            IE2Y02 = IEY02
            IE2Z02 = IEZ02
            IE2X11 = IEX11
            IE2Y11 = IEY11
            IE2Z11 = IEZ11
         END IF
         IF (IELAST .GT. LWORK2) THEN
            WRITE (LUPRI, 5000) IELAST, LWORK2
            STOP 
         END IF
      END IF
      RETURN
 1000 FORMAT(//,' <<<<<<<<<< SUBROUTINE EXCOEF >>>>>>>>>>',
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
