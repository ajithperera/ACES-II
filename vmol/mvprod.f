      SUBROUTINE MVPROD(AB,A,NA,B,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AB(NA,NB),A(NA),B(NB)
C
      DO IB=1,NB
      DO IA=1,NA
      AB(IA,IB) = A(IA)*B(IB)
      END DO
      END DO
      RETURN
C
      ENTRY PQ(AB,A,NA,B,NB,C)
      DO IB=1,NB
      D = -C*B(IB)
      DO IA=1,NA
      AB(IA,IB) = A(IA) + D
      END DO
      END DO
      RETURN
C
      ENTRY MVSUM(AB,A,NA,B,NB)
      DO IB=1,NB
      DO IA=1,NA
      AB(IA,IB) = A(IA)+B(IB)
      END DO
      END DO
      RETURN
C
      ENTRY MPPROD(AB,A,NA,B,NB,C)
      DO IB=1,NB
      D = C*B(IB)
      DO IA=1,NA
      AB(IA,IB) = A(IA)*D
      END DO
      END DO
      RETURN
C
      ENTRY MVPLUS(AB,A,NA,B,NB)
      DO IB=1,NB
      DO IA=1,NA
      AB(IA,IB) = AB(IA,IB) + A(IA)*B(IB)
      END DO
      END DO
      RETURN
C
      END
