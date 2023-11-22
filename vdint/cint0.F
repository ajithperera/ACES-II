      SUBROUTINE CINT0(NHKTA,NHKTB,KHKTA,KHKTB,ISTEPA,ISTEPB,
     *                 ISTEPU,ISTEPV,NAHGTF,NATOMC,LDIAG,IPRINT)
C
C     TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (D0 = 0.0 D00, DP5 = 0.5 D00)
      LOGICAL LDIAG
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "lwork.par"
#include "mxcent.par"
      PARAMETER (LWKRLM = LWORK3 - 8020)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      COMMON /CWORK3/ WK3LOW, SHGTF, AHGTF((MXCENT+1)*(2*MXQN+1)**3)
CSSS     &                RLMCOF(LWKRLM),WK3HGH
#include "odcs.com"
      COMMON /SDER/ SDER0 (MXAQNS),
     *              SDERX (MXAQNS), SDERY (MXAQNS), SDERZ (MXAQNS),
     *              SDERXX(MXAQNS), SDERXY(MXAQNS), SDERXZ(MXAQNS),
     *              SDERYY(MXAQNS), SDERYZ(MXAQNS), SDERZZ(MXAQNS)
      COMMON /TDER/ TDER0 (MXAQNS),
     *              TDERX (MXAQNS), TDERY (MXAQNS), TDERZ (MXAQNS),
     *              TDERXX(MXAQNS), TDERXY(MXAQNS), TDERXZ(MXAQNS),
     *              TDERYY(MXAQNS), TDERYZ(MXAQNS), TDERZZ(MXAQNS)
      COMMON /ADER/ ADER0 (MXAQNS),
     *              IA0000, IA0X00, IA0Y00, IA0Z00,
     *              IAXX00, IAXY00, IAXZ00, IAYY00,
     *              IAYZ00, IAZZ00, IA000X, IA000Y,
     *              IA000Z, IA00XX, IA00XY, IA00XZ,
     *              IA00YY, IA00YZ, IA00ZZ, IA0X0X,
     *              IA0X0Y, IA0X0Z, IA0Y0X, IA0Y0Y,
     *              IA0Y0Z, IA0Z0X, IA0Z0Y, IA0Z0Z
      COMMON /DDER/ XINT0(MXAQNS),
     *              XINTX(MXAQNS), XINTY(MXAQNS), XINTZ(MXAQNS),
     *              YINT0(MXAQNS),
     *              YINTX(MXAQNS), YINTY(MXAQNS), YINTZ(MXAQNS),
     *              ZINT0(MXAQNS),
     *              ZINTX(MXAQNS), ZINTY(MXAQNS), ZINTZ(MXAQNS),
     *              SINT0(MXAQNS)
      COMMON /RDER/ IR00, IR0X, IR0Y, IR0Z,
     *              IRXX, IRXY, IRXZ, IRYY, IRYZ, IRZZ
      COMMON /LMNS/ LVALUA(MXAQN), MVALUA(MXAQN), NVALUA(MXAQN),
     *              LVALUB(MXAQN), MVALUB(MXAQN), NVALUB(MXAQN)
      SFAC = SHGTF**3
      TFAC = -DP5*SFAC
      INT = 0
      DO 100 ICOMPA = 1,KHKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
         ISTRAT = 1 + ISTEPA*LVALA
         ISTRAU = 1 + ISTEPA*MVALA
         ISTRAV = 1 + ISTEPA*NVALA
         IF (LDIAG) THEN
            ICMPMX = MIN(ICOMPA,KHKTB)
         ELSE
            ICMPMX = KHKTB
         END IF
      DO 100 ICOMPB = 1,ICMPMX
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
         ISTRET = ISTRAT + ISTEPB*LVALB
         ISTREU = ISTRAU + ISTEPB*MVALB
         ISTREV = ISTRAV + ISTEPB*NVALB
C
C    **********************************************************
C    ***** CALCULATE OVERLAP AND KINETIC ENERGY INTEGRALS *****
C    **********************************************************
C
         X0 = ODC00X(ISTRET)
         Y0 = ODC00Y(ISTREU)
         Z0 = ODC00Z(ISTREV)
         X2 = ODC20X(ISTRET)
         Y2 = ODC20Y(ISTREU)
         Z2 = ODC20Z(ISTREV)
         INT = INT + 1
         SDER0(INT) = SDER0(INT) + SFAC*X0*Y0*Z0
         TDER0(INT) = TDER0(INT) + TFAC*(X2*Y0*Z0 + X0*Y2*Z0 + X0*Y0*Z2)
C
C     **************************************************
C     ***** CALCULATE NUCLEAR ATTRACTION INTEGRALS *****
C     **************************************************
C
         IADRAV = 1
         AINT = D0
         DO 200 IV = 0, NVALA + NVALB
            EV = ODC00Z(ISTREV + IV)
            IADRAU = IADRAV
            DO 300 IU = 0, MVALA + MVALB
               EE = ODC00Y(ISTREU + IU)*EV
               DO 400 IT = 0, LVALA + LVALB
                  EEE = ODC00X(ISTRET + IT)*EE
                  IADR00 = IADRAU + IT
                  IADD = - NAHGTF
                  DO 500 IATOM = 1,NATOMC
                     IADD = IADD + NAHGTF
                     AINT = AINT + EEE*AHGTF(IADR00 + IADD)
  500             CONTINUE
  400          CONTINUE
               IADRAU = IADRAU + ISTEPU
  300       CONTINUE
            IADRAV = IADRAV + ISTEPV
  200    CONTINUE
         ADER0(INT) = ADER0(INT) + AINT
C
C        PRINT SECTION
C
         IF (IPRINT .GE. 10) THEN
            CALL HEADER ('OUTPUT FROM CINT0',-1)
            WRITE (LUPRI,'(//,A,2I5)') ' ICOMPA, ICOMPB ', ICOMPA,ICOMPB
            WRITE (LUPRI,'(/,A,2F12.6)') ' SFAC, TFAC ', SFAC, TFAC
            WRITE (LUPRI,'(/,A,3F12.6)') ' SINT, TINT, AINT ',
     *        SFAC*X0*Y0*Z0, TFAC*(X2*Y0*Z0 + X0*Y2*Z0 + X0*Y0*Z2), AINT
         END IF
  100 CONTINUE
      RETURN
      END
