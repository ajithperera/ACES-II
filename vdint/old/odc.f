      SUBROUTINE ODC(NHKTA,NHKTB,ISTEPA,ISTEPB,FAC,EXPPI,SIGN,
     *               DIFPAX,DIFPBX,DIFPAY,DIFPBY,DIFPAZ,DIFPBZ,
     *               SIGNAX,SIGNBX,SIGNAY,SIGNBY,SIGNAZ,SIGNBZ,
     *               IPRINT)
C
C     TUH
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -------------------------------------------------------
      PARAMETER (LUCMD = 5, LUPRI = 6)
C
C     This subroutine calculates the expansion coefficients E, D,
C     and F defined by McMurchie and Davidson in J. Comp. Phys.
C     26 (1978) 218, EQS. (2.17), (2.23), AND (2.24). The recursion
C     Relations (2.20) and (2.21) are used.
C
C     The final coefficients are found in array ODC00 as illustrated
C     by the following example:
C
C     NHKTA = 2 (P) (input parameter)
C     NHKTB = 3 (D) (input parameter)
C     LA = NHKTA - 1 = 1
C     LB = NHKTB - 1 = 2
C     ISTEPB = NHKTA + NHKTB = 5 (input parameter)
C     ISTEPA = NHKTB*ISTEPB = 15 (input parameter)
C
C     E(IA,IB,IT)
C     IA = 0, LA
C     IB = 0, LB
C     IT = 0, IA + IB
C
C     E(0,0,0)    0.0       0.0       0.0       0.0
C     E(0,1,0)  E(0,1,1)    0.0       0.0       0.0
C     E(0,2,0)  E(0,2,1)  E(0,2,2)    0.0       0.0
C     E(1,0,0)  E(1,0,1)    0.0       0.0       0.0
C     E(1,1,0)  E(1,1,1)  E(1,1,2)    0.0       0.0
C     E(1,2,0)  E(1,2,1)  E(1,2,2)  E(1,2,3)    0.0
C
C     Address of element E(0,0,0) = FAC is 1.
C
C     Dimension requirements on ODC00(-I:J) are
C
C     I = MAX(NHKTB*(NHKTA + NHKTB))
C     J = MAX(NHKTA*NHKTB*(NHKTA + NHKTB))
C
      PARAMETER (ZERO = 0.00 D00, HALF = 0.50 D00)
      COMMON /ODCS/ ODC10X(192), ODC20X(192), ODC30X(192), ODC40X(192),
     *              ODC00X(432), ODC10Y(192), ODC20Y(192), ODC30Y(192),
     *              ODC40Y(192), ODC00Y(432), ODC10Z(192), ODC20Z(192),
     *              ODC30Z(192), ODC40Z(192), ODC00Z(432)
      LA = NHKTA - 1
      LB = NHKTB - 1
      DO 100 I = 0,NHKTA*ISTEPA
         ODC00X(I) = ZERO
         ODC00Y(I) = ZERO
         ODC00Z(I) = ZERO
  100 CONTINUE
      ISTART = 1
      LAST0 = ISTART - ISTEPA
      LASTM = LAST0 - 1
      LASTP = LAST0 + 1
      ODC00X(LAST0) = ZERO
      ODC00Y(LAST0) = ZERO
      ODC00Z(LAST0) = ZERO
      ODC00X(LASTM) = ZERO
      ODC00Y(LASTM) = ZERO
      ODC00Z(LASTM) = ZERO
      FAC = SIGN*FAC
      ODC00X(LASTP) = SIGNAX*FAC
      ODC00Y(LASTP) = SIGNAY*FAC
      ODC00Z(LASTP) = SIGNAZ*FAC
      EXPPIH = SIGN*HALF*EXPPI
      DO 200 IA = 0,LA
         DO 300 IT = 0,IA
            T1 = SIGN*FLOAT(IT + 1)
            IADR  = ISTART + IT
            LAST0 = IADR - ISTEPA
            LASTM = LAST0 - 1
            LASTP = LAST0 + 1
            ODC00X(IADR) = SIGNAX*(EXPPIH*ODC00X(LASTM)
     *                           + DIFPAX*ODC00X(LAST0)
     *                           + T1*ODC00X(LASTP))
            ODC00Y(IADR) = SIGNAY*(EXPPIH*ODC00Y(LASTM)
     *                           + DIFPAY*ODC00Y(LAST0)
     *                           + T1*ODC00Y(LASTP))
            ODC00Z(IADR) = SIGNAZ*(EXPPIH*ODC00Z(LASTM)
     *                           + DIFPAZ*ODC00Z(LAST0)
     *                           + T1*ODC00Z(LASTP))
  300    CONTINUE
         ISTART = ISTART + ISTEPB
         DO 400 IB = 1,LB
            DO 500 IT = 0,IA + IB
               T1 = SIGN*FLOAT(IT + 1)
               IADR  = ISTART + IT
               LAST0 = IADR - ISTEPB
               LASTM = LAST0 - 1
               LASTP = LAST0 + 1
               ODC00X(IADR) = SIGNBX*(EXPPIH*ODC00X(LASTM)
     *                              + DIFPBX*ODC00X(LAST0)
     *                              + T1*ODC00X(LASTP))
               ODC00Y(IADR) = SIGNBY*(EXPPIH*ODC00Y(LASTM)
     *                              + DIFPBY*ODC00Y(LAST0)
     *                              + T1*ODC00Y(LASTP))
               ODC00Z(IADR) = SIGNBZ*(EXPPIH*ODC00Z(LASTM)
     *                              + DIFPBZ*ODC00Z(LAST0)
     *                              + T1*ODC00Z(LASTP))
  500       CONTINUE
            ISTART = ISTART + ISTEPB
  400    CONTINUE
  200 CONTINUE
      IF (IPRINT .GE. 20) THEN
         CALL TITLER('Output from ODC','*',103)
         NHKTAB = NHKTA*NHKTB
         CALL AROUND('ODC00X')
         CALL OUTPUT(ODC00X,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC00Y')
         CALL OUTPUT(ODC00Y,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
         CALL AROUND('ODC00Z')
         CALL OUTPUT(ODC00Z,1,NHKTAB,1,ISTEPB,NHKTAB,ISTEPB,1,LUPRI)
      END IF
      RETURN
      END
