
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

      SUBROUTINE ODC(NHKTA,  NHKTB,
     &               ISTEPA, ISTEPB,
     &               FAC,
     &               EXPPI,  SIGN,
     &               DIFPAX, DIFPBX,
     &               DIFPAY, DIFPBY,
     &               DIFPAZ, DIFPBZ,
     &               SIGNAX, SIGNBX,
     &               SIGNAY, SIGNBY,
     &               SIGNAZ, SIGNBZ,
     &               IPRINT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision dtmp(3)
      PARAMETER (LUCMD = 5, LUPRI = 6)
      PARAMETER (MXQN=8, MXAQN=MXQN*(MXQN+1)/2, MXAQNS=MXAQN*MXAQN)
#include "odcs.com" 
c ----------------------------------------------------------------------

      ITMP = NHKTA*ISTEPA
      DO I = 0, ITMP
         ODC00X(I) = 0.0d0
      END DO
      DO I = 0, ITMP
         ODC00Y(I) = 0.0d0
      END DO
      DO I = 0, ITMP
         ODC00Z(I) = 0.0d0
      END DO

      FAC = SIGN*FAC
      IADR = -ISTEPA
      ODC00X(IADR    ) = 0.0d0
      ODC00X(IADR + 1) = 0.0d0
      ODC00X(IADR + 2) = SIGNAX*FAC
      ODC00Y(IADR    ) = 0.0d0
      ODC00Y(IADR + 1) = 0.0d0
      ODC00Y(IADR + 2) = SIGNAY*FAC
      ODC00Z(IADR    ) = 0.0d0
      ODC00Z(IADR + 1) = 0.0d0
      ODC00Z(IADR + 2) = SIGNAZ*FAC

      dtmp(1) = SIGN*0.5d0*EXPPI

      LA = NHKTA - 1
      LB = NHKTB - 1

      ISTART = 1
      DO IA = 0, LA
         dtmp(2) = DIFPAX
         DO IT = 0, IA
            dtmp(3) = SIGN*FLOAT(IT + 1)
            IADR  = ISTART + IT
            LASTM = IADR - ( ISTEPA + 1 )
            ODC00X(IADR) =   ODC00X(LASTM    ) * dtmp(1)
     &                     + ODC00X(LASTM + 1) * dtmp(2)
     &                     + ODC00X(LASTM + 2) * dtmp(3)
            ODC00X(IADR) = ODC00X(IADR) * SIGNAX
         END DO
         ISTART = ISTART + ISTEPB
         dtmp(2) = DIFPBX
         DO IB = 1, LB
            ITMP = IA + IB
            DO IT = 0, ITMP
               dtmp(3) = SIGN*FLOAT(IT + 1)
               IADR  = ISTART + IT
               LASTM = IADR - ( ISTEPB + 1 )
               ODC00X(IADR) =   dtmp(1) * ODC00X(LASTM    )
     &                        + dtmp(2) * ODC00X(LASTM + 1)
     &                        + dtmp(3) * ODC00X(LASTM + 2)
               ODC00X(IADR) = ODC00X(IADR) * SIGNBX
            END DO
            ISTART = ISTART + ISTEPB
         END DO
      END DO

      ISTART = 1
      DO IA = 0, LA
         dtmp(2) = DIFPAY
         DO IT = 0, IA
            dtmp(3) = SIGN*FLOAT(IT + 1)
            IADR  = ISTART + IT
            LASTM = IADR - ( ISTEPA + 1 )
            ODC00Y(IADR) =   ODC00Y(LASTM    ) * dtmp(1)
     &                     + ODC00Y(LASTM + 1) * dtmp(2)
     &                     + ODC00Y(LASTM + 2) * dtmp(3)
            ODC00Y(IADR) = ODC00Y(IADR) * SIGNAY
         END DO
         ISTART = ISTART + ISTEPB
         dtmp(2) = DIFPBY
         DO IB = 1, LB
            ITMP = IA + IB
            DO IT = 0, ITMP
               dtmp(3) = SIGN*FLOAT(IT + 1)
               IADR  = ISTART + IT
               LASTM = IADR - ( ISTEPB + 1 )
               ODC00Y(IADR) =   dtmp(1) * ODC00Y(LASTM    )
     &                        + dtmp(2) * ODC00Y(LASTM + 1)
     &                        + dtmp(3) * ODC00Y(LASTM + 2)
               ODC00Y(IADR) = ODC00Y(IADR) * SIGNBY
            END DO
            ISTART = ISTART + ISTEPB
         END DO
      END DO

      ISTART = 1
      DO IA = 0, LA
         dtmp(2) = DIFPAZ
         DO IT = 0, IA
            dtmp(3) = SIGN*FLOAT(IT + 1)
            IADR  = ISTART + IT
            LASTM = IADR - ( ISTEPA + 1 )
            ODC00Z(IADR) =   ODC00Z(LASTM    ) * dtmp(1)
     &                     + ODC00Z(LASTM + 1) * dtmp(2)
     &                     + ODC00Z(LASTM + 2) * dtmp(3)
            ODC00Z(IADR) = ODC00Z(IADR) * SIGNAZ
         END DO
         ISTART = ISTART + ISTEPB
         dtmp(2) = DIFPBZ
         DO IB = 1, LB
            ITMP = IA + IB
            DO IT = 0, ITMP
               dtmp(3) = SIGN*FLOAT(IT + 1)
               IADR  = ISTART + IT
               LASTM = IADR - ( ISTEPB + 1 )
               ODC00Z(IADR) =   dtmp(1) * ODC00Z(LASTM    )
     &                        + dtmp(2) * ODC00Z(LASTM + 1)
     &                        + dtmp(3) * ODC00Z(LASTM + 2)
               ODC00Z(IADR) = ODC00Z(IADR) * SIGNBZ
            END DO
            ISTART = ISTART + ISTEPB
         END DO
      END DO
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

