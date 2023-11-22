      SUBROUTINE TODC(NHKTA,NHKTB,ISTEPA,ISTEPB,EXPA)
C
C     TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FAC21(8), FAC20(8)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)

#include "odcs.com" 
C
C     The DATA statements below are tabulations of the
C     the following functions:
C
C     FAC21(N) = 2*N + 1
C     FAC20(N) = N*(N - 1)
C
      DATA FAC21 /1.D0, 3.D0, 5.D0, 7.D0, 9.D0, 11.D0, 13.D0, 15.D0/
     *     FAC20 /0.D0, 0.D0, 2.D0, 6.D0, 12.D0, 20.D0, 30.D0, 42.D0/
      TEXPA1 = EXPA + EXPA
      TEXPA2 = TEXPA1*TEXPA1
      IADR0 = 1
      IADRP1 = IADR0  + ISTEPA
      IADRP2 = IADRP1 + ISTEPA
      IADRM1 = IADR0  - ISTEPA
      IADRM2 = IADRM1 - ISTEPA
      DO 100 IA = 1,NHKTA
         FAC21A = TEXPA1*FAC21(IA)
         FAC20A = FAC20(IA)
         DO 200 IB = 1,NHKTB
            ODC20X(IADR0) = TEXPA2*ODC00X(IADRP2)
     *                    - FAC21A*ODC00X(IADR0)
     *                    + FAC20A*ODC00X(IADRM2)
            ODC20Y(IADR0) = TEXPA2*ODC00Y(IADRP2)
     *                    - FAC21A*ODC00Y(IADR0)
     *                    + FAC20A*ODC00Y(IADRM2)
            ODC20Z(IADR0) = TEXPA2*ODC00Z(IADRP2)
     *                    - FAC21A*ODC00Z(IADR0)
     *                    + FAC20A*ODC00Z(IADRM2)
            IADR0  = IADR0  + ISTEPB
            IADRP1 = IADRP1 + ISTEPB
            IADRP2 = IADRP2 + ISTEPB
            IADRM1 = IADRM1 + ISTEPB
            IADRM2 = IADRM2 + ISTEPB
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
