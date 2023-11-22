      SUBROUTINE DODC(NHKTA,NHKTB,ISTEPA,ISTEPB,EXPA,IPRINT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
      DIMENSION FAC10(8), FAC21(8), FAC20(8), FAC32(8),
     *          FAC31(8), FAC30(8), FAC43(8), FAC42(8),
     *          FAC41(8), FAC40(8)
#include "odcs.com"
C
C     The DATA statements below are tabulations of the
C     the following functions:
C
C     FAC10(N) = N
C     FAC21(N) = 2*N + 1
C     FAC20(N) = N*(N - 1)
C     FAC32(N) = 3*N + 3
C     FAC31(N) = 3*N*N
C     FAC30(N) = N*(N - 1)*(N - 2)
C     FAC43(N) = 4*N + 6
C     FAC42(N) = 6*N*(N + 1) + 3
C     FAC41(N) = N*(N - 1)*(4*N - 2)
C     FAC40(N) = N*(N - 1)*(N - 2)*(N - 3)
C
      DATA FAC10 /0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0/
     &     FAC21 /1.D0,3.D0,5.D0,7.D0,9.D0,11.D0,13.D0,15.D0/
     &     FAC20 /0.D0,0.D0,2.D0,6.D0,12.D0,20.D0,30.D0,42.D0/
     &     FAC32 /3.D0,6.D0,9.D0,12.D0,15.D0,18.D0,21.D0,24.D0/
     &     FAC31 /0.D0,3.D0,12.D0,27.D0,48.D0,75.D0,108.D0,147.D0/
     &     FAC30 /0.D0,0.D0,0.D0,6.D0,24.D0,60.D0,120.D0,210.D0/
     &     FAC43 /6.D0,10.D0,14.D0,18.D0,22.D0,26.D0,30.D0,34.D0/
     &     FAC42 /3.D0,15.D0,39.D0,75.D0,123.D0,183.D0,255.D0,339.D0/
     &     FAC41 /0.D0,0.D0,12.D0,60.D0,168.D0,360.D0,660.D0,1092.D0/
     &     FAC40 /0.D0,0.D0,0.D0,0.D0,24.D0,120.D0,360.D0,840.D0/

      call zero(odc10x,1280) 
      call zero(odc10y,1280) 
      call zero(odc10z,1280) 
      call zero(odc20x,1280) 
      call zero(odc20y,1280) 
      call zero(odc20z,1280) 
      TEXPA1 = EXPA + EXPA
      TEXPA2 = TEXPA1*TEXPA1
      TEXPA3 = TEXPA2*TEXPA1
      TEXPA4 = TEXPA3*TEXPA1
      IADR0 = 1
      IADRP1 = IADR0  + ISTEPA
      IADRP2 = IADRP1 + ISTEPA
      IADRP3 = IADRP2 + ISTEPA
      IADRP4 = IADRP3 + ISTEPA
      IADRM1 = IADR0  - ISTEPA
      IADRM2 = IADRM1 - ISTEPA
      IADRM3 = IADRM2 - ISTEPA
      IADRM4 = IADRM3 - ISTEPA
      DO 100 IA = 1,NHKTA
         FAC10A = FAC10(IA)
         FAC21A = TEXPA1*FAC21(IA)
         FAC20A = FAC20(IA)
         FAC32A = TEXPA2*FAC32(IA)
         FAC31A = TEXPA1*FAC31(IA)
         FAC30A = FAC30(IA)
         FAC43A = TEXPA3*FAC43(IA)
         FAC42A = TEXPA2*FAC42(IA)
         FAC41A = TEXPA1*FAC41(IA)
         FAC40A = FAC40(IA)
         DO 200 IB = 1,NHKTB
            DO 300 I = 0, IA + IB
               ODC10X(IADR0 + I) = TEXPA1*ODC00X(IADRP1 + I)
     *                           - FAC10A*ODC00X(IADRM1 + I)
               ODC10Y(IADR0 + I) = TEXPA1*ODC00Y(IADRP1 + I)
     *                           - FAC10A*ODC00Y(IADRM1 + I)
               ODC10Z(IADR0 + I) = TEXPA1*ODC00Z(IADRP1 + I)
     *                           - FAC10A*ODC00Z(IADRM1 + I)
               ODC20X(IADR0 + I) = TEXPA2*ODC00X(IADRP2 + I)
     *                           - FAC21A*ODC00X(IADR0  + I)
     *                           + FAC20A*ODC00X(IADRM2 + I)
               ODC20Y(IADR0 + I) = TEXPA2*ODC00Y(IADRP2 + I)
     *                           - FAC21A*ODC00Y(IADR0  + I)
     *                           + FAC20A*ODC00Y(IADRM2 + I)
               ODC20Z(IADR0 + I) = TEXPA2*ODC00Z(IADRP2 + I)
     *                           - FAC21A*ODC00Z(IADR0  + I)
     *                           + FAC20A*ODC00Z(IADRM2 + I)
  300       CONTINUE
            ODC30X(IADR0) = TEXPA3*ODC00X(IADRP3)
     *                    - FAC32A*ODC00X(IADRP1)
     *                    + FAC31A*ODC00X(IADRM1)
     *                    - FAC30A*ODC00X(IADRM3)
            ODC30Y(IADR0) = TEXPA3*ODC00Y(IADRP3)
     *                    - FAC32A*ODC00Y(IADRP1)
     *                    + FAC31A*ODC00Y(IADRM1)
     *                    - FAC30A*ODC00Y(IADRM3)
            ODC30Z(IADR0) = TEXPA3*ODC00Z(IADRP3)
     *                    - FAC32A*ODC00Z(IADRP1)
     *                    + FAC31A*ODC00Z(IADRM1)
     *                    - FAC30A*ODC00Z(IADRM3)
            ODC40X(IADR0) = TEXPA4*ODC00X(IADRP4)
     *                    - FAC43A*ODC00X(IADRP2)
     *                    + FAC42A*ODC00X(IADR0)
     *                    - FAC41A*ODC00X(IADRM2)
     *                    + FAC40A*ODC00X(IADRM4)
            ODC40Y(IADR0) = TEXPA4*ODC00Y(IADRP4)
     *                    - FAC43A*ODC00Y(IADRP2)
     *                    + FAC42A*ODC00Y(IADR0)
     *                    - FAC41A*ODC00Y(IADRM2)
     *                    + FAC40A*ODC00Y(IADRM4)
            ODC40Z(IADR0) = TEXPA4*ODC00Z(IADRP4)
     *                    - FAC43A*ODC00Z(IADRP2)
     *                    + FAC42A*ODC00Z(IADR0)
     *                    - FAC41A*ODC00Z(IADRM2)
     *                    + FAC40A*ODC00Z(IADRM4)
            IADR0 = IADR0 + ISTEPB
            IADRP1 = IADRP1 + ISTEPB
            IADRP2 = IADRP2 + ISTEPB
            IADRP3 = IADRP3 + ISTEPB
            IADRP4 = IADRP4 + ISTEPB
            IADRM1 = IADRM1 + ISTEPB
            IADRM2 = IADRM2 + ISTEPB
            IADRM3 = IADRM3 + ISTEPB
            IADRM4 = IADRM4 + ISTEPB
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from DODC','*',103)
         NHKTAB = NHKTA*NHKTB
         CALL AROUND('ODC10X')
         CALL OUTPUT(ODC10X,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC20X')
         CALL OUTPUT(ODC20X,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC30X')
         CALL OUTPUT(ODC30X,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC40X')
         CALL OUTPUT(ODC40X,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC10Y')
         CALL OUTPUT(ODC10Y,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC20Y')
         CALL OUTPUT(ODC20Y,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC30Y')
         CALL OUTPUT(ODC30Y,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC40Y')
         CALL OUTPUT(ODC40Y,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC10Z')
         CALL OUTPUT(ODC10Z,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC20Z')
         CALL OUTPUT(ODC20Z,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC30Z')
         CALL OUTPUT(ODC30Z,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC40Z')
         CALL OUTPUT(ODC40Z,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
      END IF
      RETURN
      END
