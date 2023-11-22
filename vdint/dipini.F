      SUBROUTINE DIPINI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "baslims.par"
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2,
     *           MXAQNS=MXAQN*MXAQN*MXCONT*MXCONT)
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
      CALL ZERO(XINT0, (13*MXAQNS) )
      RETURN
      END
