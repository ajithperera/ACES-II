      SUBROUTINE MVPROD(AB,A,NA,B,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AB(NA,NB),A(NA),B(NB)
C
      IF(NA.LT.NB) THEN
      DO 1 IA=1,NA
      DO 1 IB=1,NB
    1 AB(IA,IB) = A(IA)*B(IB)
C
      ELSE
      DO 2 IB=1,NB
      DO 2 IA=1,NA
    2 AB(IA,IB) = A(IA)*B(IB)
      ENDIF
C
      RETURN
      ENTRY PQ (AB,A,NA,B,NB,C)
      IF(NA.LT.NB) THEN
C
      DO 3 IA=1,NA
      DO 3 IB=1,NB
    3 AB(IA,IB) =A(IA) - C*B(IB)
C
      ELSE
      DO 4 IB=1,NB
      DO 4 IA=1,NA
    4 AB(IA,IB) = A(IA) - C*B(IB)
C
      ENDIF
      RETURN
      ENTRY MVSUM (AB,A,NA,B,NB)
C
      IF(NA.LT.NB) THEN
      DO 5 IA=1,NA
      DO 5 IB=1,NB
    5 AB(IA,IB) = A(IA)+B(IB)
C
      ELSE
      DO 6 IB=1,NB
      DO 6 IA=1,NA
    6 AB(IA,IB) = A(IA)+B(IB)
C
      ENDIF
      RETURN
      ENTRY MPPROD (AB,A,NA,B,NB,C)
      IF(NA.LT.NB) THEN
      DO 7 IA=1,NA
      DO 7 IB=1,NB
    7 AB(IA,IB) = A(IA)*B(IB)*C
C
      ELSE
      DO 8 IB=1,NB
      DO 8 IA=1,NA
    8 AB(IA,IB) = A(IA)*B(IB)*C
      ENDIF
      RETURN
C
      ENTRY MVPLUS (AB,A,NA,B,NB)
C
      IF(NA.LT.NB) THEN
      DO 9 IA=1,NA
      DO 9 IB=1,NB
    9 AB(IA,IB) = AB(IA,IB) + A(IA)*B(IB)
C
      ELSE
      DO 10 IB=1,NB
      DO 10 IA=1,NA
   10 AB(IA,IB) = AB(IA,IB) + A(IA)*B(IB)
C
      ENDIF
      RETURN
      END
