      

      SUBROUTINE GIADRV(NHKT1,NHKT2,NSET1,NSET2,ISTEP1,ISTEP2,
     &                   NUC1,NUC2,NUC12,NORB1,NORB2,
     &                   NUCO1,NUCO2,NRCO1,NRCO2,JSTR1,JSTR2,
     &                   SIGN1X,SIGN1Y,SIGN1Z,SIGN2X,SIGN2Y,SIGN2Z,
     &                   COR1X,COR1Y,COR1Z,COR2X,COR2Y,COR2Z,
     &                   SAMEX,SAMEY,SAMEZ,I0X,I0Y,I0Z,TPRIAB,
     &                   ONECEN,DO1,DO2,BIGVEC,SEGMEN,DTEST,
     &                   ITYPE,THRESH,MAXDER,IPRINT)
C
CEND
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
#include "baslims.par"
      PARAMETER (ONE = 1.00 D00)
      LOGICAL SAMEX, SAMEY, SAMEZ, ONECEN, DO1, DO2, TPRIAB,
     &        BIGVEC, SEGMEN, DTEST
      DIMENSION NUCO1(MXAOVC), NUCO2(MXAOVC),
     &          NRCO1(MXCONT), NRCO2(MXCONT),
     &          JSTR1(MXAOVC), JSTR2(MXAOVC)
C
      I0X = 0
      I0Y = 0
      I0Z = 0
C
C OVERLAP DISTRIBUTION VECTORS 
C
      CALL GIAVEC(NUC1,NUC2,NORB1,NORB2,NSET1,NSET2,NUCO1,NUCO2,
     &             NRCO1,NRCO2,JSTR1,JSTR2,NUC12,TPRIAB,
     &             SIGN1X,SIGN1Y,SIGN1Z,
     &             SIGN2X,SIGN2Y,SIGN2Z,
     &             THRESH,ITYPE,IPRINT)
      IF (NUC12 .EQ. 0) RETURN
C
C OVERLAP DISTRIBUTION COEFFICIENTS 
C
      CALL GIACOE(BIGVEC,ONECEN,SEGMEN,NHKT1,NHKT2,ISTEP1,ISTEP2,NUC12,
     &            I0X,I0Y,I0Z,DIFX,DIFY,DIFZ,SIGN1X,SIGN1Y,SIGN1Z,
     &            SIGN2X,SIGN2Y,SIGN2Z,THRESH,ITYPE,MAXDER,
     &            DO1,DO2,DTEST,IPRINT)
      RETURN
      END
